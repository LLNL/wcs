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
#include "utils/rngen.hpp"
#include "utils/trace_ssa.hpp"
#include "utils/timer.hpp"
#include "sim_methods/ssa_nrm.hpp"
#include "sim_methods/ssa_direct.hpp"


#define OPTIONS "dg:hi:o:s:t:m:r:"
static const struct option longopts[] = {
    {"diag",     no_argument,        0, 'd'},
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
    time_interval(0.0)
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
    "    -d, --diag\n"
    "            Specify whether to enable tracing for a posteriori diagnosis.\n"
    "\n"
    "    -g, --graphviz\n"
    "            Specify the name of the file to export the reaction\n"
    "            network into in the GraphViz format.\n"
    "\n"
    "    -h, --help\n"
    "            Display this usage information\n"
    "\n"
    "    -o, --outfile\n"
    "            Specify the output file name for tracing/sampling.\n"
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
    "            Specify the SSA method: 0 = direct, 1 = next reaction (default).\n"
    "\n"
    "    -r, --record\n"
    "            Specify whether to enable sampling at a time/step interval.\n"
    "            (e.g. t5 for every 5 simulation seconds or i5 for every 5 steps).\n"
    "\n";
  exit(code);
}


int main(int argc, char** argv)
{
  int rc = 0;
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
    rc = -1;
  }

  wcs::Sim_Method* ssa = nullptr;

  if (cfg.method == 0) {
    ssa = new wcs::SSA_Direct;
    std::cerr << "Direct SSA method." << std::endl;
  } else if (cfg.method == 1) {
    std::cerr << "Next Reaction SSA method." << std::endl;
    ssa = new wcs::SSA_NRM;
  } else {
    std::cerr << "Unknown SSA method (" << cfg.method << ')' << std::endl;
    return -1;
  }
  if (cfg.tracing) {
    ssa->set_tracing();
    std::cerr << "Enable tracing" << std::endl;
  } else if (cfg.sampling) {
    if (cfg.iter_interval > 0u) {
      ssa->set_sampling(cfg.iter_interval);
      std::cerr << "Enable sampling at " << cfg.iter_interval 
                << " steps interval" << std::endl;
    } else {
      ssa->set_sampling(cfg.time_interval);
      std::cerr << "Enable sampling at " << cfg.time_interval 
                << " secs interval" << std::endl;
    }
  }
  ssa->init(rnet_ptr, cfg.max_iter, cfg.max_time, cfg.seed);
  double t_start = wcs::get_time();
  ssa->run();
  std::cout << "Wall clock time to run simulation: "
            << wcs::get_time() - t_start << " (sec)" << std::endl;

  if (cfg.tracing) {
    if (!cfg.outfile.empty()) {
      ssa->trace().write(cfg.outfile);
    } else {
      ssa->trace().write(std::cout);
    }
  } else if (cfg.sampling) {
    if (!cfg.outfile.empty()) {
      ssa->samples().write(cfg.outfile);
    } else {
      ssa->samples().write(std::cout);
    }
  } else {
    std::cout << "Species   : " << rnet.show_species_labels("") << std::endl;
    std::cout << "FinalState: " << rnet.show_species_counts() << std::endl;
  }

  delete ssa;

  return rc;
}
