/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include <getopt.h>
#include <fstream>
#include <cmath>
#include <limits>
#include <vector>
#include <set>
#include <string>
#include <boost/filesystem.hpp>
#include "reaction_network/network.hpp"
#include "utils/write_graphviz.hpp"
#include "utils/timer.hpp"
#include "sim_methods/ssa_nrm.hpp"
#include "sim_methods/ssa_direct.hpp"
#include "sim_methods/ssa_sod.hpp"

#ifdef WCS_HAS_VTUNE
__itt_domain* vtune_domain_sim = __itt_domain_create("Simulate");
__itt_string_handle* vtune_handle_sim = __itt_string_handle_create("simulate");
#endif // WCS_HAS_VTUNE

#define OPTIONS "df:g:hi:o:s:t:m:r:"
static const struct option longopts[] = {
    {"diag",     no_argument,        0, 'd'},
    {"frag_sz",  required_argument,  0, 'f'},
    {"graphviz", required_argument,  0, 'g'},
    {"help",     no_argument,        0, 'h'},
    {"iter",     required_argument,  0, 'i'},
    {"outfile",  required_argument,  0, 'o'},
    {"seed",     required_argument,  0, 's'},
    {"time",     required_argument,  0, 't'},
    {"method",   required_argument,  0, 'm'},
    {"record",   required_argument,  0, 'r'},
    { 0, 0, 0, 0 },
};

struct Config {
  Config()
  : seed(0u), max_iter(10u),
    max_time(wcs::Network::get_etime_ulimit()),
    method(1),
    tracing(false),
    sampling(false),
    iter_interval(0u),
    time_interval(0.0),
    frag_size(0),
    is_frag_size_set(false)
  {}

  void getopt(int& argc, char** &argv);
  void print_usage(const std::string exec, int code);

  unsigned seed;
  wcs::sim_iter_t max_iter;
  wcs::sim_time_t max_time;
  int method;
  bool tracing;
  bool sampling;
  wcs::sim_iter_t iter_interval;
  wcs::sim_time_t time_interval;
  unsigned frag_size;
  bool is_frag_size_set;

  std::string infile;
  std::string outfile;
  std::string gvizfile;
};

void Config::getopt(int& argc, char** &argv)
{
  int c;
  bool is_iter_set = false;
  bool is_time_set = false;

  while ((c = getopt_long(argc, argv, OPTIONS, longopts, NULL)) != -1) {
    switch (c) {
      case 'd': /* --diag */
        tracing = true;
        sampling = false;
        break;
      case 'f': /* --frag_sz */
        frag_size = static_cast<unsigned>(atoi(optarg));
        is_frag_size_set = true;
        break;
      case 'g': /* --graphviz */
        gvizfile = std::string(optarg);
        break;
      case 'h': /* --help */
        print_usage(argv[0], 0);
        break;
      case 'i': /* --iter */
        max_iter = static_cast<wcs::sim_iter_t>(atoi(optarg));
        is_iter_set = true;
        break;
      case 'o': /* --outfile */
        outfile = std::string(optarg);
        break;
      case 's': /* --seed */
        seed = static_cast<unsigned>(atoi(optarg));
        break;
      case 't': /* --time */
        max_time = static_cast<wcs::sim_time_t>(std::stod(optarg));
        is_time_set = true;
        break;
      case 'm': /* --method */
        method = static_cast<int>(atoi(optarg));
        break;
      case 'r': /* --record */
        sampling = true;
        tracing = false;
        {
          if (optarg[0] == 'i') {
            iter_interval = static_cast<wcs::sim_iter_t>(atoi(&optarg[1]));
          } else if (optarg[0] == 't') {
            time_interval = static_cast<wcs::sim_time_t>(atof(&optarg[1]));
          } else {
            sampling = false;
            std::cerr << "Unknown sampling interval: " << std::string(optarg)
                      << std::endl;
          }
        }
        break;
      default:
        print_usage(argv[0], 1);
        break;
    }
  }

  if (optind != (argc - 1)) {
    print_usage (argv[0], 1);
  }

  infile = argv[optind];

  if (!is_iter_set && is_time_set) {
    max_iter = std::numeric_limits<decltype(max_iter)>::max();
  }
  if ((sampling || tracing) && !is_frag_size_set) {
    frag_size = wcs::default_frag_size;
  }
}

void Config::print_usage(const std::string exec, int code)
{
  std::cerr <<
    "Usage: " << exec << " <filename>.graphml\n"
    "    Run the stochastic simulation algorithm (SSA) on a reaction\n"
    "    network graph (<filename>.graphml) for the given number of\n"
    "    iteration and with the given random number seed.\n"
    "    Specifically, the next reaction method is used. Optionally,\n"
    "    write the reaction network into a GraphViz-formatted file.\n"
    "\n"
    "    OPTIONS:\n"
    "    -h, --help\n"
    "            Display this usage information\n"
    "\n"
    "    -s, --seed\n"
    "            Specify the seed for random number generator. Without this,\n"
    "            it will use a value dependent on the current system clock.\n"
    "\n"
    "    -t, --time\n"
    "            Specify the upper limit of simulation time to run.\n"
    "\n"
    "    -i, --iter\n"
    "            Specify the maximum number of reaction iterations to run.\n"
    "\n"
    "    -m, --method\n"
    "            Specify the SSA method: 0 = Optimized direct method "
    " (feat. dependency graph)\n"
    "                                    1 = next reaction (default).\n"
    "                                    2 = Sorted optimized direct method."
    " (feat. propensity sorting)\n"
    "\n"
    "    -g, --graphviz\n"
    "            Specify the name of the file to export the reaction\n"
    "            network into in the GraphViz format.\n"
    "\n"
    "    -d, --diag\n"
    "            Specify whether to enable tracing for a posteriori diagnosis.\n"
    "\n"
    "    -r, --record\n"
    "            Specify whether to enable sampling at a time/step interval.\n"
    "            (e.g. t5 for every 5 simulation seconds or i5 for every 5 steps).\n"
    "\n"
    "    -o, --outfile\n"
    "            Specify the output file name for tracing/sampling.\n"
    "\n"
    "    -f, --frag_sz\n"
    "            Specify how many records per temporary output file fragment \n"
    "            in tracing/sampling.\n"
    "\n";
  exit(code);
}


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

  if (cfg.method == 0) {
    ssa = new wcs::SSA_Direct;
    std::cerr << "Direct SSA method." << std::endl;
  } else if (cfg.method == 1) {
    std::cerr << "Next Reaction SSA method." << std::endl;
    ssa = new wcs::SSA_NRM;
  } else if (cfg.method == 2) {
    std::cerr << "Sorted optimized direct SSA method." << std::endl;
    ssa = new wcs::SSA_SOD;
  } else {
    std::cerr << "Unknown SSA method (" << cfg.method << ')' << std::endl;
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
  ssa->init(rnet_ptr, cfg.max_iter, cfg.max_time, cfg.seed);

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
