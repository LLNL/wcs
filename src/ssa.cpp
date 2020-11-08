/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include <string>
#include <iostream>
#include "ssa-cfg.hpp"
#include "utils/write_graphviz.hpp"
#include "utils/timer.hpp"
#include "reaction_network/network.hpp"
#include "sim_methods/ssa_nrm.hpp"
#include "sim_methods/ssa_direct.hpp"
#include "sim_methods/ssa_sod.hpp"

#ifdef WCS_HAS_VTUNE
__itt_domain* vtune_domain_sim = __itt_domain_create("Simulate");
__itt_string_handle* vtune_handle_sim = __itt_string_handle_create("simulate");
#endif // WCS_HAS_VTUNE


int main(int argc, char** argv)
{
 #ifdef WCS_HAS_VTUNE
  __itt_pause();
 #endif // WCS_HAS_VTUNE
  int rc = EXIT_SUCCESS;
  Config cfg;
  cfg.getopt(argc, argv);

  std::shared_ptr<wcs::Network> rnet_ptr = std::make_shared<wcs::Network>();
  wcs::Network& rnet = *rnet_ptr;
  rnet.load(cfg.infile);
  rnet.init();
  const wcs::Network::graph_t& g = rnet.graph();

  if (!cfg.gvizfile.empty() &&
      !wcs::write_graphviz(cfg.gvizfile, g))
  {
    std::cerr << "Failed to write " << cfg.gvizfile << std::endl;
    rc = EXIT_FAILURE;
  }

  wcs::Sim_Method* ssa = nullptr;

  try {
    if (cfg.method == 0) {
      ssa = new wcs::SSA_Direct(rnet_ptr);
      std::cerr << "Direct SSA method." << std::endl;
    } else if (cfg.method == 1) {
      std::cerr << "Next Reaction SSA method." << std::endl;
      ssa = new wcs::SSA_NRM(rnet_ptr);
    } else if (cfg.method == 2) {
      std::cerr << "Sorted optimized direct SSA method." << std::endl;
      ssa = new wcs::SSA_SOD(rnet_ptr);
    } else {
      std::cerr << "Unknown SSA method (" << cfg.method << ')' << std::endl;
      return EXIT_FAILURE;
    }
  } catch (const std::exception& e) {
    std::cerr << "Fail to setup SSA method." << std::endl;
    return EXIT_FAILURE;
  }

  if (cfg.tracing) {
    ssa->set_tracing<wcs::TraceSSA>(cfg.outfile, cfg.frag_size);
    std::cerr << "Enable tracing" << std::endl;
  } else if (cfg.sampling) {
    if (cfg.iter_interval > 0u) {
      ssa->set_sampling<wcs::SamplesSSA>(cfg.iter_interval, cfg.outfile, cfg.frag_size);
      std::cerr << "Enable sampling at " << cfg.iter_interval
                << " steps interval" << std::endl;
    } else {
      ssa->set_sampling<wcs::SamplesSSA>(cfg.time_interval, cfg.outfile, cfg.frag_size);
      std::cerr << "Enable sampling at " << cfg.time_interval
                << " secs interval" << std::endl;
    }
  }
  ssa->init(cfg.max_iter, cfg.max_time, cfg.seed);

 #ifdef WCS_HAS_VTUNE
  __itt_resume();
  __itt_task_begin(vtune_domain_sim, __itt_null, __itt_null, vtune_handle_sim);
 #endif // WCS_HAS_VTUNE

  double t_start = wcs::get_time();
  ssa->run();
  std::cout << "Wall clock time to run simulation: "
            << wcs::get_time() - t_start << " (sec)" << std::endl;

 #ifdef WCS_HAS_VTUNE
  __itt_task_end(vtune_domain_sim);
  __itt_pause();
 #endif // WCS_HAS_VTUNE

  if (cfg.tracing || cfg.sampling) {
    ssa->finalize_recording();
  } else {
    std::cout << "Species   : " << rnet.show_species_labels("") << std::endl;
    std::cout << "FinalState: " << rnet.show_species_counts() << std::endl;
  }

  delete ssa;

  return rc;
}
