/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#if defined(WCS_HAS_ROSS)

#include <string>
#include <cstring> // memset
#include <iostream>
#include "utils/write_graphviz.hpp"
#include "utils/timer.hpp"
#include "utils/to_string.hpp"
#include "reaction_network/network.hpp"
#include "des.hpp"
#include "wcs-ross-bf.hpp"

using revent_t = wcs::Sim_State_Change::revent_t;
WCS_Global_State gState;

int main(int argc, char **argv)
{
  tw_opt_add(app_opt);
  tw_init(&argc, &argv);

  wcs::SSA_Params& cfg = gState.m_cfg;
  cfg.getopt(argc, argv);
  tw_define_lps(nlp_per_pe, sizeof(WCS_Message));

  for(unsigned int i = 0; i < g_tw_nlp; i++) {
    tw_lp_settype(i, &wcs_LPs[0]);
  }

  if( g_tw_mynode == 0 )
  {
    std::cout << "=========================================" << std::endl;
    std::cout << "WCS ROSS Configuration.............." << std::endl;
    std::cout << "   run_id:\t" + std::string(run_id) << std::endl;
    std::cout << "   nlp_per_pe:\t" << g_tw_nlp << std::endl;
    std::cout << "   g_tw_ts_end:\t" << g_tw_ts_end << std::endl;;
    std::cout << "   gvt-interval:\t" << g_tw_gvt_interval << std::endl;;
    std::cout << "   extramem:\t" << g_tw_events_per_pe_extra << std::endl;
    std::cout << "   ......................................" << std::endl;
    std::cout << "   Num nodes:\t" << tw_nnodes() << std::endl;
    std::cout << "   Message size:\t" << sizeof(WCS_Message) << std::endl;
    std::cout << "========================================="
              << std::endl << std::endl;
  }

  tw_run();
  tw_end();

  return EXIT_SUCCESS;
}


void wcs_init(WCS_State *s, tw_lp *lp)
{
  wcs::SSA_Params& cfg = gState.m_cfg;

  std::shared_ptr<wcs::Network> rnet_ptr = std::make_shared<wcs::Network>();
  wcs::Network& rnet = *rnet_ptr;
  rnet.load(cfg.m_infile);
  rnet.init();
  const wcs::Network::graph_t& g = rnet.graph();

  if (!cfg.m_gvizfile.empty() &&
      !wcs::write_graphviz(cfg.m_gvizfile, g))
  {
    WCS_THROW("Failed to write " + cfg.m_gvizfile);
    return;
  }

  std::unique_ptr<wcs::SSA_NRM> ssa;

  try {
    if (cfg.m_method == 1) {
      std::cerr << "Next Reaction SSA method." << std::endl;
      ssa = std::make_unique<wcs::SSA_NRM>(rnet_ptr);
    } else {
      WCS_THROW("Unsupported SSA method (" + std::to_string(cfg.m_method) + ')');
      return;
    }
  } catch (const std::exception& e) {
    WCS_THROW("Fail to setup SSA method.");
    return;
  }

  if (cfg.m_tracing) {
    ssa->set_tracing<wcs::TraceSSA>(cfg.get_outfile(), cfg.m_frag_size);
    std::cerr << "Enable tracing" << std::endl;
  } else if (cfg.m_sampling) {
    if (cfg.m_iter_interval > 0u) {
      ssa->set_sampling<wcs::SamplesSSA>(cfg.m_iter_interval,
                                         cfg.get_outfile(),
                                         cfg.m_frag_size);
      std::cerr << "Enable sampling at " << cfg.m_iter_interval
                << " steps interval" << std::endl;
    } else {
      ssa->set_sampling<wcs::SamplesSSA>(cfg.m_time_interval,
                                         cfg.get_outfile(),
                                         cfg.m_frag_size);
      std::cerr << "Enable sampling at " << cfg.m_time_interval
                << " secs interval" << std::endl;
    }
  }
  ssa->init(cfg.m_max_iter, cfg.m_max_time, cfg.m_seed);

  const size_t lp_idx = gState.m_LP_states.size();
  ssa->m_lp_idx = lp_idx;
  s->m_lp_idx = lp_idx;
  gState.m_LP_states.emplace_back(std::move(ssa), rnet_ptr);
  gState.m_LP_states[lp_idx].m_t_start = wcs::get_time();
}


void wcs_prerun(WCS_State *s, tw_lp *lp)
{
  const WCS_LP_State& lp_state = gState.m_LP_states.at(s->m_lp_idx);

  wcs::Sim_Method::revent_t reaction_1st;

  if (lp_state.m_ssa_ptr->schedule(reaction_1st) == wcs::Sim_Method::Success)
  {
    tw_event* next_evt = tw_event_new(lp->gid, reaction_1st.first, lp);
    auto* next_msg = reinterpret_cast<WCS_Message*>(tw_event_data(next_evt));
    next_msg->reaction = lp_state.m_net_ptr->reaction_d2i(reaction_1st.second);
    tw_event_send(next_evt);
  }
}


void wcs_event(WCS_State *s, tw_bf *bf, WCS_Message *msg, tw_lp *lp)
{
  memset(static_cast<void*>(bf), 0, sizeof(tw_bf));
  const WCS_LP_State& lp_state = gState.m_LP_states.at(s->m_lp_idx);

  wcs::Sim_Method::revent_t firing
    = std::make_pair(tw_now(lp), lp_state.m_net_ptr->reaction_i2d(msg->reaction));

  if (lp_state.m_ssa_ptr->forward(firing))
  {
    WCS_BF_(bf, WCS_BF_FWD) = 1u;
    wcs::SSA_NRM::priority_t new_firing;
    if (lp_state.m_ssa_ptr->schedule(new_firing) == wcs::Sim_Method::Success)
    {
      WCS_BF_(bf, WCS_BF_SCHED) = 1u;
      new_firing.first -= tw_now(lp);
      tw_event* next_evt = tw_event_new(lp->gid, new_firing.first, lp);
      auto* next_msg = reinterpret_cast<WCS_Message*>(tw_event_data(next_evt));
      next_msg->reaction = lp_state.m_net_ptr->reaction_d2i(new_firing.second);
      tw_event_send(next_evt);
    }
  }
}


void wcs_event_reverse(WCS_State *s, tw_bf *bf, WCS_Message *msg, tw_lp *lp)
{
  if (!WCS_BF_(bf, WCS_BF_FWD)) {
    return;
  }
  const WCS_LP_State& lp_state = gState.m_LP_states.at(s->m_lp_idx);

  wcs::Sim_Method::revent_t firing
    = std::make_pair(tw_now(lp), lp_state.m_net_ptr->reaction_i2d(msg->reaction));

  lp_state.m_ssa_ptr->backward(firing);
}


void wcs_event_commit(WCS_State *s, tw_bf *bf, WCS_Message *msg, tw_lp *lp)
{
  if (!WCS_BF_(bf, WCS_BF_FWD)) {
    return;
  }

  const WCS_LP_State& lp_state = gState.m_LP_states.at(s->m_lp_idx);

  lp_state.m_ssa_ptr->commit_des();
}


void wcs_final(WCS_State *s, tw_lp *lp)
{
  wcs::SSA_Params& cfg = gState.m_cfg;
  const WCS_LP_State& lp_state = gState.m_LP_states.at(s->m_lp_idx);

  std::cout << "Wall clock time to run simulation: "
            << wcs::get_time() - lp_state.m_t_start << " (sec)" << std::endl;

  if (cfg.m_tracing || cfg.m_sampling) {
    lp_state.m_ssa_ptr->finalize_recording();
  } else {
    std::cout << "Species   : "
              << lp_state.m_net_ptr->show_species_labels("") << std::endl;
    std::cout << "FinalState: "
              << lp_state.m_net_ptr->show_species_counts() << std::endl;
  }
}


tw_peid wcs_map(tw_lpid gid)
{
  return (tw_peid) gid / g_tw_nlp;
}

#endif // defined(WCS_HAS_ROSS)
