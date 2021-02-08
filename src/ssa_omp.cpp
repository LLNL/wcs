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
#include "utils/timer.hpp"
#include "utils/write_graphviz.hpp"
#include "params/wcs_params.hpp"
#include "params/ssa_params.hpp"
#include "proto/wcs_params.hpp"
#include "partition/metis_partition.hpp"
#include "partition/partition_info.hpp"
#include "sim_methods/ssa_nrm.hpp"

/// Shared state
struct WCS_Shared_State {
  wcs::sim_time_t m_max_time;
  wcs::sim_iter_t m_max_iter;
  int m_nparts;
};

WCS_Shared_State shared_state;

/// logical process state
struct WCS_LP_State
{
  /// Pointer to an ssa object
  std::unique_ptr<wcs::SSA_NRM> m_ssa_ptr;
  /// Pointer to the reaction network
  std::shared_ptr<wcs::Network> m_net_ptr;
  /// Simulation start time (wall clock)
  double m_t_start;

  WCS_LP_State(std::unique_ptr<wcs::SSA_NRM>&& ssa_ptr,
               std::shared_ptr<wcs::Network>& net_ptr);

  WCS_LP_State() = default;
  WCS_LP_State(WCS_LP_State&& other) = default;
  WCS_LP_State& operator=(WCS_LP_State&& other) = default;

};

WCS_LP_State::WCS_LP_State(std::unique_ptr<wcs::SSA_NRM>&& ssa_ptr,
                           std::shared_ptr<wcs::Network>& net_ptr)
: m_ssa_ptr(std::move(ssa_ptr)), m_net_ptr(net_ptr), m_t_start(0.0)
{}

namespace {
/// Thread-local variable, file-visible only.
#ifdef __ICC
  WCS_LP_State lp_state;
  #pragma omp threadprivate(lp_state)
#else
  // Defined like this to work around a GCC problem with threadprivate objects:
  // https://stackoverflow.com/questions/23552077/how-to-define-a-object-or-struct-as-threadprivate-in-openmp/
  extern WCS_LP_State lp_state;
  #pragma omp threadprivate(lp_state)
  WCS_LP_State lp_state;
#endif
}

/// Type to exchange a reaction betwen different subdomains on different LPs
using nrm_evt_t = std::pair<wcs::sim_time_t, wcs::v_idx_t>;

/// Result of partitioning
using partition_idx_t = std::vector<idx_t>;

/// Partition the given reaction network using Metis
bool initial_partition(const wcs::Metis_Params& mp,
                       partition_idx_t& parts,
                       idx_t& objval,
                       bool write_subgraph = false);

bool setup_partition(const std::string& input_model,
                     const wcs::Metis_Params& mp_in,
                     partition_idx_t& parts);

void wcs_init(const wcs::SSA_Params& cfg, const partition_idx_t& parts);

std::pair<wcs::sim_iter_t, wcs::sim_time_t> wcs_run();


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

  partition_idx_t parts; ///< Partition assignment result
  ok = setup_partition(cmd.m_input_model, mp, parts);
  if (!ok) return EXIT_FAILURE;

  google::protobuf::ShutdownProtobufLibrary();
  sp.m_infile = cmd.m_input_model;
  if (sp.m_method == 1) {
    std::cerr << "Next Reaction SSA method." << std::endl;
  } else {
    WCS_THROW("Ony NRM is supported");
    return EXIT_FAILURE;
  }

 // ......................... SIMULATION BEGINS ................................
 #if defined(_OPENMP) && defined(WCS_OMP_RUN_PARTITION)
   wcs_init(sp, parts);

 #ifdef WCS_HAS_VTUNE
  __itt_resume();
  __itt_task_begin(vtune_domain_sim, __itt_null, __itt_null, vtune_handle_sim);
 #endif // WCS_HAS_VTUNE

  double t_start = wcs::get_time();
  wcs_run();
  std::cout << "Wall clock time to run simulation: "
            << wcs::get_time() - t_start << " (sec)" << std::endl;

 #ifdef WCS_HAS_VTUNE
  __itt_task_end(vtune_domain_sim);
  __itt_pause();
 #endif // WCS_HAS_VTUNE

  #pragma omp master
  {
    if (sp.m_tracing || sp.m_sampling) {
      auto& ssa = *(lp_state.m_ssa_ptr);
      ssa.finalize_recording();
    } else {
      const auto& net = *(lp_state.m_net_ptr);
      std::cout << "Species   : " << net.show_species_labels("") << std::endl;
      std::cout << "FinalState: " << net.show_species_counts() << std::endl;
    }
  }
 #endif // defined(_OPENMP) && defined(WCS_OMP_RUN_PARTITION)

  return EXIT_SUCCESS;
}


//-----------------------------------------------------------------------------
#if defined(_OPENMP) && defined(WCS_OMP_RUN_PARTITION)
//-----------------------------------------------------------------------------
void wcs_init(const wcs::SSA_Params& cfg, const partition_idx_t& parts)
{
  shared_state.m_max_time = cfg.m_max_time;
  shared_state.m_max_iter = cfg.m_max_iter;

  #pragma omp parallel num_threads(shared_state.m_nparts)
  {
    auto& ssa_ptr = lp_state.m_ssa_ptr;
    auto& net_ptr = lp_state.m_net_ptr;

    net_ptr = std::make_shared<wcs::Network>();
    wcs::Network& net = *net_ptr;

    net.load(cfg.m_infile);
    net.init();
    net.set_partition(parts, omp_get_thread_num());

    ssa_ptr = std::make_unique<wcs::SSA_NRM>(net_ptr);
    wcs::SSA_NRM& ssa = *ssa_ptr;

    if (cfg.m_tracing) {
      ssa.set_tracing<wcs::TraceSSA>(cfg.m_outfile, cfg.m_frag_size);
      std::cerr << "Enable tracing" << std::endl;
    } else if (cfg.m_sampling) {
      if (cfg.m_iter_interval > 0u) {
        ssa.set_sampling<wcs::SamplesSSA>(cfg.m_iter_interval,
                                          cfg.m_outfile,
                                          cfg.m_frag_size);
        std::cerr << "Enable sampling at " << cfg.m_iter_interval
                  << " steps interval" << std::endl;
      } else {
        ssa.set_sampling<wcs::SamplesSSA>(cfg.m_time_interval,
                                          cfg.m_outfile,
                                          cfg.m_frag_size);
        std::cerr << "Enable sampling at " << cfg.m_time_interval
                  << " secs interval" << std::endl;
      }
    }
    ssa.init(cfg.m_max_iter, cfg.m_max_time, cfg.m_seed);

    lp_state.m_t_start = wcs::get_time();
  }
}


wcs::Sim_Method::result_t schedule(nrm_evt_t& evt_earliest)
{
  constexpr nrm_evt_t sevt_undef {std::numeric_limits<wcs::sim_time_t>::max(),
                                  std::numeric_limits<wcs::v_idx_t>::max()};

  evt_earliest = sevt_undef;

  #pragma omp declare reduction(earliest :\
    nrm_evt_t :\
    omp_out = \
      (((omp_in.first < omp_out.first) || \
        ((omp_in.first == omp_out.first) && (omp_in.second < omp_out.second)))? \
       omp_in : omp_out)) \
    initializer(omp_priv = sevt_undef)

  #pragma omp parallel num_threads(shared_state.m_nparts) reduction(earliest:evt_earliest)
  {
    const auto& ssa = *(lp_state.m_ssa_ptr);
    const auto& net = *(lp_state.m_net_ptr);

    if (BOOST_UNLIKELY(ssa.is_empty())) {
      evt_earliest = sevt_undef;
    } else {
      wcs::Sim_Method::revent_t re = ssa.choose_reaction();

      if (BOOST_UNLIKELY(re.first > shared_state.m_max_time)) {
        evt_earliest = sevt_undef;
      } else {
        evt_earliest = nrm_evt_t{re.first, net.reaction_d2i(re.second)};
      }
    }
  }

  if (BOOST_UNLIKELY(evt_earliest.first > shared_state.m_max_time)) {
    return wcs::Sim_Method::Inactive;
  }
  return wcs::Sim_Method::Success;
}


bool forward(const nrm_evt_t& evt_earliest)
{
  bool ok[omp_get_max_threads()] = {false};
  #pragma omp parallel num_threads(shared_state.m_nparts)
  {
    auto& ssa = *(lp_state.m_ssa_ptr);
    ok[omp_get_thread_num()] = ssa.advance_time_and_iter(evt_earliest.first);
  }
  if (!ok[0]) return false;

  #pragma omp parallel num_threads(shared_state.m_nparts)
  {
    auto& ssa = *(lp_state.m_ssa_ptr);
    const auto& net = *(lp_state.m_net_ptr);

    wcs::Sim_Method::revent_t firing =
      std::make_pair (evt_earliest.first,
                      net.reaction_i2d(evt_earliest.second));

    wcs::Sim_State_Change digest(firing);

    // Depdending on whether the firing reaction is local or not,
    // the processing of the reaction is different.
    const bool local = (net.graph()[firing.second].get_partition() ==
                        omp_get_thread_num());
    if (local) {
      // Execute the reaction, updating species counts
      ssa.fire_reaction(digest);
      // update the propensities and times of those reactions fired and affected
      ssa.update_reactions(firing, digest.m_reactions_affected, digest.m_reaction_times);
    } else {
      // Only returns the affected reactions that are local
      ssa.fire_reaction(digest);
      // Don't update the reaction fired which is not local.
      ssa.update_reactions(firing.first, digest.m_reactions_affected, digest.m_reaction_times);
    }

    #pragma omp master
    {
      ssa.record(firing.second);
    }
  }
  return true;
}


std::pair<wcs::sim_iter_t, wcs::sim_time_t> wcs_run()
{
  nrm_evt_t next_reaction;

  if (schedule(next_reaction) != wcs::Sim_Method::Success) {
    WCS_THROW("Not able to schedule any reaction event!");
  }

  while (BOOST_LIKELY(forward(next_reaction))) {
    if (BOOST_UNLIKELY(schedule(next_reaction) != wcs::Sim_Method::Success)) {
      break;
    }
  }

  std::pair<wcs::sim_iter_t, wcs::sim_time_t> ret;
  #pragma omp master
  {
    const auto& ssa = *(lp_state.m_ssa_ptr);
    ret = std::make_pair(ssa.get_sim_iter(), ssa.get_sim_time());
  }
  return ret;
}

//-----------------------------------------------------------------------------
#endif // if defined(_OPENMP) && defined(WCS_OMP_RUN_PARTITION)
//-----------------------------------------------------------------------------


bool initial_partition(const wcs::Metis_Params& mp,
                       partition_idx_t& parts,
                       idx_t& objval,
                       bool write_subgraph)
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

    if (mp.m_verbose && !(mp.m_rnet->my_reaction_list()).empty()) {
      if (write_subgraph) {
        const auto gpart_name
          = wcs::append_to_stem(mp.m_outfile, "-" + std::to_string(i));
        if (!wcs::write_graphviz(gpart_name, mp.m_rnet->graph(), i)) {
          std::cerr << "Failed to write " << gpart_name << std::endl;
          continue;
        }
      }
      std::cout << "(" + std::to_string((mp.m_rnet->my_reaction_list()).size())
        + " reactions, " + std::to_string((mp.m_rnet->my_species_list()).size())
        + " species)" << std::endl;
    }
  }
  return true;
}

bool setup_partition(const std::string& input_model,
                     const wcs::Metis_Params& mp_in,
                     partition_idx_t& parts)
{
  std::shared_ptr<wcs::Network> rnet_ptr = std::make_shared<wcs::Network>();
  wcs::Network& rnet = *rnet_ptr;

  rnet.load(input_model);
  rnet.init();
  if (rnet.get_num_reactions() == 0u) {
    std::cerr << "No reaction available in network!" << std::endl;
    return false;
  }

  wcs::Metis_Params mp = mp_in;
  mp.m_rnet = rnet_ptr;

  idx_t objval; /// Total comm volume or edge-cut of the solution
  bool ok = initial_partition(mp, parts, objval);
  if (!ok) return false;

  if (mp.m_verbose) {
    wcs::Partition_Info pinfo(rnet_ptr);
    pinfo.scan(mp.m_verbose);
    pinfo.report();
  }

  shared_state.m_nparts = mp.m_nparts;
  return true;
}

#else  // defined(WCS_HAS_METIS)
#error This code requires METIS
#endif // defined(WCS_HAS_METIS)
