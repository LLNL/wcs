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
#include "utils/omp_diagnostics.hpp"
#include "params/wcs_params.hpp"
#include "params/ssa_params.hpp"
#include "proto/wcs_params.hpp"
#include "partition/metis_partition.hpp"
#include "partition/partition_info.hpp"
#include "sim_methods/ssa_nrm.hpp"

#ifdef WCS_HAS_VTUNE
__itt_domain* vtune_domain_sim = __itt_domain_create("Simulate");
__itt_string_handle* vtune_handle_sim = __itt_string_handle_create("simulate");
#endif // WCS_HAS_VTUNE

/// Shared state
struct WCS_Shared_State {
  wcs::sim_time_t m_max_time;
  wcs::sim_iter_t m_max_iter;
  int m_nparts;
  int m_num_inner_threads;

  void set_num_partitions(int np);
};

void WCS_Shared_State::set_num_partitions(int np)
{
  m_nparts = np;
  if (m_nparts <= 0) {
    m_nparts = 1;
    WCS_THROW("Invalid number of partitions: " + std::to_string(np));
  }
 #if defined(_OPENMP)
  m_num_inner_threads = omp_get_max_threads() / m_nparts;
  if (m_num_inner_threads <= 0) {
    std::string msg = "Number of partitions (" + std::to_string(m_nparts);
    msg += ") must be less than or equal to the number of available threads (";
    msg += std::to_string(omp_get_max_threads());
    WCS_THROW(msg);
  }
 #else
  m_num_inner_threads = 1;
 #endif // defined(_OPENMP)
}

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

#if defined(_OPENMP) && defined(WCS_OMP_RUN_PARTITION)
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
#endif // defined(_OPENMP) && defined(WCS_OMP_RUN_PARTITION)

/// Type to exchange a reaction betwen different subdomains on different LPs
using nrm_evt_t = std::pair<wcs::sim_time_t, wcs::v_idx_t>;

/// Result of partitioning
using partition_idx_t = std::vector<idx_t>;

/// Partition the given reaction network using Metis
bool setup_partition(const std::string& input_model,
                     const wcs::Metis_Params& mp_in,
                     partition_idx_t& parts,
                     const bool write_subgraph = false);

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
  if (sp.m_method != 1) {
    WCS_THROW("Ony NRM is supported");
    return EXIT_FAILURE;
  }

 // ......................... SIMULATION BEGINS ................................
 #if defined(_OPENMP) && defined(WCS_OMP_RUN_PARTITION)
  wcs_init(sp, parts);

  std::cout << "Initialization complete. Simulation begins ..." << std::endl;

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
      std::string ofile = sp.get_outfile();
      if (ofile.empty()) {
        ofile = wcs::get_default_outname_from_model(cmd.m_input_model);
      }
      std::ofstream ofs(ofile);
      ofs << "Species   : " << net.show_species_labels("") << std::endl;
      ofs << "FinalState: " << net.show_species_counts() << std::endl;
    }
  }
 #elif defined(_OPENMP)
  std::cout << "This mode of parallelization does not require network "
            << "partitioning. Use `ssa` instead." << std::endl;
 #endif // defined(_OPENMP) && defined(WCS_OMP_RUN_PARTITION)

  return EXIT_SUCCESS;
}


//-----------------------------------------------------------------------------
#if defined(_OPENMP) && defined(WCS_OMP_RUN_PARTITION)
//-----------------------------------------------------------------------------
void wcs_init(const wcs::SSA_Params& cfg, const partition_idx_t& parts)
{
#if defined(WCS_OMP_REACTION_UPDATES) || \
    defined(WCS_OMP_REACTION_REACTANTS) || \
    defined(WCS_OMP_REACTION_PRODUCTS)
  //omp_set_nested(1);
  omp_set_max_active_levels(2);
#endif
  omp_set_dynamic(0);
  omp_set_schedule(omp_sched_static, 0);
  omp_set_num_threads(shared_state.m_nparts);

  shared_state.m_max_time = cfg.m_max_time;
  shared_state.m_max_iter = cfg.m_max_iter;

 #if OMP_DEBUG
  std::vector<wcs::my_omp_affinity> omp_aff(shared_state.m_nparts);
 #endif // OMP_DEBUG

  #pragma omp parallel num_threads(shared_state.m_nparts)
  {
   #if OMP_DEBUG
    omp_aff[omp_get_thread_num()].get();
   #endif // OMP_DEBUG

    auto& ssa_ptr = lp_state.m_ssa_ptr;
    auto& net_ptr = lp_state.m_net_ptr;

    net_ptr = std::make_shared<wcs::Network>();
    wcs::Network& net = *net_ptr;

    net.load(cfg.m_infile, true);
    net.init();
    net.set_partition(parts, omp_get_thread_num());

    ssa_ptr = std::make_unique<wcs::SSA_NRM>(net_ptr);
    wcs::SSA_NRM& ssa = *ssa_ptr;
    ssa.set_num_threads(shared_state.m_num_inner_threads);

    #pragma omp master
    {
      if (cfg.m_tracing) {
        ssa.set_tracing<wcs::TraceSSA>(cfg.get_outfile(), cfg.m_frag_size);
      } else if (cfg.m_sampling) {
        if (cfg.m_iter_interval > 0u) {
          ssa.set_sampling<wcs::SamplesSSA>(cfg.m_iter_interval,
                                            cfg.get_outfile(),
                                            cfg.m_frag_size);
        } else {
          ssa.set_sampling<wcs::SamplesSSA>(cfg.m_time_interval,
                                            cfg.get_outfile(),
                                            cfg.m_frag_size);
        }
      }
    }
    ssa.init(cfg.m_max_iter, cfg.m_max_time, cfg.m_seed);
    ssa.m_lp_idx = omp_get_thread_num();

    lp_state.m_t_start = wcs::get_time();
  }

 #if OMP_DEBUG
  for (const auto& oaff: omp_aff) {
    oaff.print();
  }
 #endif // OMP_DEBUG

  if (cfg.m_tracing) {
      std::cerr << "Enable tracing" << std::endl;
  } else if (cfg.m_sampling) {
    if (cfg.m_iter_interval > 0u) {
      std::cerr << "Enable sampling at " << cfg.m_iter_interval
                << " steps interval" << std::endl;
    } else {
      std::cerr << "Enable sampling at " << cfg.m_time_interval
                << " secs interval" << std::endl;
    }
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


bool setup_partition(const std::string& input_model,
                     const wcs::Metis_Params& mp_in,
                     partition_idx_t& parts,
                     const bool write_subgraph)
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

  wcs::Metis_Partition partitioner(mp);
  partitioner.prepare();

  //const auto& map_idx2vd = partitioner.get_map_from_idx_to_desc();
  if (mp.m_verbose) {
    partitioner.print_params();
    partitioner.print_metis_graph(std::cout);
    partitioner.print_adjacency(std::cout);
  }

  idx_t objval; /// Total comm volume or edge-cut of the solution
  bool ret = partitioner.run(parts, objval);
  if (!ret) return false;
  shared_state.set_num_partitions(mp.m_nparts);

  if (mp.m_verbose) {
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
    wcs::Partition_Info pinfo(rnet_ptr);
    pinfo.scan(mp.m_verbose);
    pinfo.report();
  }
  return true;
}

#else  // defined(WCS_HAS_METIS)
#error This code requires METIS
#endif // defined(WCS_HAS_METIS)
