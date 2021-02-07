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

#if !defined(WCS_HAS_PROTOBUF)
#error WCS requires protocol buffer
#endif

// This code requires METIS
#if defined(WCS_HAS_METIS)
#include <string>
#include <iostream>
#include <vector>
#include <functional>
#include "utils/file.hpp"
#include "utils/write_graphviz.hpp"
#include "params/wcs_params.hpp"
#include "params/ssa_params.hpp"
#include "proto/wcs_params.hpp"
#include "partition/metis_partition.hpp"
#include "partition/partition_info.hpp"
#include "sim_methods/ssa_nrm.hpp"

/// Partition the given reaction network using Metis
bool initial_partition(const wcs::Metis_Params& mp,
                       std::vector<idx_t>& parts,
                       idx_t& objval)
{
  wcs::Metis_Partition partitioner(mp);
  partitioner.prepare();
  //const auto& map_idx2vd = partitioner.get_map_from_idx_to_desc();
  if (mp.m_verbose) {
    partitioner.print_params();
    partitioner.print_metis_graph(std::cout);
    partitioner.print_adjacency(std::cout);
  }

  bool ret = partitioner.run(parts, objval);
  if (!ret) return false;

  const auto& map_idx2desc = partitioner.get_map_from_idx_to_desc();

  for (wcs::partition_id_t i = 0; i < mp.m_nparts; ++i) {
    mp.m_rnet->set_partition(map_idx2desc, parts, i);

    if (!(mp.m_rnet->my_reaction_list()).empty()) {
      const auto gpart_name
        = wcs::append_to_stem(mp.m_outfile, "-" + std::to_string(i));

      wcs::SSA_NRM nrm(mp.m_rnet);
      nrm.init(1, 0.0, mp.get_seed());
      const auto r = nrm.choose_reaction();
      const auto rname = ((mp.m_rnet->graph())[r.second]).get_label();
      std::cout << "First reaction from partition " << i << " is "
                << rname << " at time "  << r.first << std::endl;

      if (!wcs::write_graphviz(gpart_name, mp.m_rnet->graph(), i)) {
        std::cerr << "Failed to write " << gpart_name << std::endl;
        continue;
      }
      std::cout << "(" + std::to_string((mp.m_rnet->my_reaction_list()).size())
        + " reactions, " + std::to_string((mp.m_rnet->my_species_list()).size())
        + " species)" << std::endl;
    }
  }
  return true;
}


int main(int argc, char** argv)
{
  wcs::cmd_line_opts cmd;
  bool ok = cmd.parse_cmd_line(argc, argv);
  if (!ok) return EXIT_FAILURE;
  if (!cmd.m_is_set) return EXIT_SUCCESS;

  cmd.show();

  wcs::SSA_Params sp;
  wcs::Metis_Params mp;

  if (!cmd.m_all_setup.empty()) {
    wcs::read_proto_params(cmd.m_all_setup, sp, mp);
  } else {
    if (!cmd.m_sim_setup.empty()) {
      wcs::read_proto_params(cmd.m_sim_setup, sp);
    }
    if (!cmd.m_part_setup.empty()) {
      wcs::read_proto_params(cmd.m_part_setup, mp);
    }
  }
  google::protobuf::ShutdownProtobufLibrary();

  std::shared_ptr<wcs::Network> rnet_ptr = std::make_shared<wcs::Network>();
  wcs::Network& rnet = *rnet_ptr;

  rnet.load(cmd.m_input_model);
  rnet.init();

  mp.m_rnet = rnet_ptr;
  sp.m_infile = cmd.m_input_model;

  std::vector<idx_t> parts; ///< Partition assignment result
  idx_t objval; /// Total comm volume or edge-cut of the solution
  ok = initial_partition(mp, parts, objval);
  if (!ok) return EXIT_FAILURE;

  wcs::Partition_Info pinfo(rnet_ptr);
  pinfo.scan(mp.m_verbose);
  pinfo.report();

  return EXIT_SUCCESS;
}
#else  // defined(WCS_HAS_METIS)
#error This code requires METIS
#endif // defined(WCS_HAS_METIS)

#if 0
void wcs_init(SSA_)
{
  wcs::SSA_Params& cfg = gState.m_cfg;

  std::shared_ptr<wcs::Network> rnet_ptr = std::make_shared<wcs::Network>();
  wcs::Network& rnet = *rnet_ptr;
  rnet.load(cfg.infile);
  rnet.init();
  const wcs::Network::graph_t& g = rnet.graph();

  if (!cfg.gvizfile.empty() &&
      !wcs::write_graphviz(cfg.gvizfile, g))
  {
    WCS_THROW("Failed to write " + cfg.gvizfile);
    return;
  }

  std::unique_ptr<wcs::SSA_NRM> ssa;

  try {
    if (cfg.method == 1) {
      std::cerr << "Next Reaction SSA method." << std::endl;
      ssa = std::make_unique<wcs::SSA_NRM>(rnet_ptr);
    } else {
      WCS_THROW("Unsupported SSA method (" + std::to_string(cfg.method) + ')');
      return;
    }
  } catch (const std::exception& e) {
    WCS_THROW("Fail to setup SSA method.");
    return;
  }

  if (cfg.tracing) {
    ssa->set_tracing<wcs::TraceSSA>(cfg.outfile, cfg.frag_size);
    std::cerr << "Enable tracing" << std::endl;
  } else if (cfg.sampling) {
    if (cfg.iter_interval > 0u) {
      ssa->set_sampling<wcs::SamplesSSA>(cfg.iter_interval,
                                         cfg.outfile,
                                         cfg.frag_size);
      std::cerr << "Enable sampling at " << cfg.iter_interval
                << " steps interval" << std::endl;
    } else {
      ssa->set_sampling<wcs::SamplesSSA>(cfg.time_interval,
                                         cfg.outfile,
                                         cfg.frag_size);
      std::cerr << "Enable sampling at " << cfg.time_interval
                << " secs interval" << std::endl;
    }
  }
  ssa->init(cfg.max_iter, cfg.max_time, cfg.seed);

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
#endif
