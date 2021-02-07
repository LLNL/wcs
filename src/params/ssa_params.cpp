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
#include <limits>
#include <string>
#include <iostream>
#include <cstdlib>
#include "reaction_network/network.hpp"
#include "params/ssa_params.hpp"


namespace wcs {

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

SSA_Params::SSA_Params()
: seed(0u), max_iter(10u),
  max_time(wcs::max_sim_time),
  method(1),
  tracing(false),
  sampling(false),
  iter_interval(0u),
  time_interval(0.0),
  frag_size(0),
  is_frag_size_set(false),
  is_iter_set(false),
  is_time_set(false)
{}

void SSA_Params::getopt(int& argc, char** &argv)
{
  int c;
  is_iter_set = false;
  is_time_set = false;

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

void SSA_Params::print_usage(const std::string exec, int code)
{
  std::cerr <<
    "Usage: " << exec << " <filename>.graphml|xml\n"
    "    Run the stochastic simulation algorithm (SSA) on a reaction\n"
    "    network graph (<filename>.graphml|xml) for the given number of\n"
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

} // end of namespace wcs
