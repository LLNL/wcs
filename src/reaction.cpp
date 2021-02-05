/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include <iostream>
#include <getopt.h>
#include <boost/filesystem.hpp>
#include "utils/write_graphviz.hpp"
#include "reaction_network/network.hpp"
#include <fstream>

#ifdef WCS_HAS_VTUNE
__itt_domain* vtune_domain_react = __itt_domain_create("Reactions");
__itt_string_handle* vtune_handle_react = __itt_string_handle_create("react");
#endif // WCS_HAS_VTUNE


#define OPTIONS "hi:o:"
static const struct option longopts[] = {
    {"help",    no_argument,  0, 'h'},
    {"num_iter", required_argument, 0, 'i'},
    {"outfile", required_argument,  0, 'o'},
    { 0, 0, 0, 0 },
};

void print_usage(const std::string exec, int code)
{
  std::cerr <<
    "Usage: " << exec << " <filename>.graphml\n"
    "    Load a reaction network graph (<filename>.graphml) or\n"
    "    a Systems Biology Markup Language (SBML) file (<filename>.xml)\n"
    "    into a graph with the flat vertex property stucture.\n"
    "    Then, convert it to the one with polymorphic vertices.\n"
    "    Finally, write to a GraphViz format file (<filename>.dot)\n"
    "    that has the same name as the input file except the \n"
    "    extention '.dot' unless --outfile is given.\n"
    "\n"
    "    OPTIONS:\n"
    "    -h, --help\n"
    "            Display this usage information\n"
    "\n"
    "    -i, --num_iter\n"
    "            Specify the number of iterations to compute all"
    "            the reaction rates for measuring time\n"
    "\n"
    "    -o, --outfile\n"
    "            Specify the output file name\n"
    "\n";
  exit(code);
}


void count_active_reactions(const wcs::Network& rnet)
{
  size_t num_active = 0ul;
  size_t num_total = rnet.reaction_list().size();

  for (const auto r : rnet.reaction_list()) {
    num_active += static_cast<size_t>(rnet.check_reaction(r));
  }

  std::cout << "Num active reactions: "
            << num_active << "/" << num_total << std::endl;
}

int main(int argc, char** argv)
{
 #ifdef WCS_HAS_VTUNE
  __itt_pause();
 #endif // WCS_HAS_VTUNE
 #ifdef WCS_HAS_HPCTOOLKIT
  hpctoolkit_sampling_stop();
 #endif // WCS_HAS_HPCTOOLKIT

  int c;
  int rc = 0;
  std::string outfile;
  unsigned int num_iter = 0u;

  while ((c = getopt_long(argc, argv, OPTIONS, longopts, NULL)) != -1) {
    switch (c) {
      case 'h': /* --help */
        print_usage(argv[0], 0);
        break;
      case 'i': /* --num_iter */
        num_iter = static_cast<unsigned>(atoi(optarg));
        break;
      case 'o': /* --outfile */
        outfile = std::string(optarg);
        break;
      default:
        print_usage(argv[0], 1);
        break;
    }
  }

  if (optind != (argc - 1)) {
    print_usage (argv[0], 1);
  }

  std::string fn(argv[optind]);

  wcs::Network rnet;

  rnet.load(fn);
  rnet.init();

  count_active_reactions(rnet);

 #ifndef WCS_PERF_PROF
  const wcs::Network::graph_t& g = rnet.graph();

  if (outfile.empty()) {
    boost::filesystem::path path = fn;
    std::string base = path.stem().string();
    outfile = base + "_poly.dot";;;
  }

  if (!wcs::write_graphviz(outfile, g)) {
    std::cerr << "Failed to write " << outfile << std::endl;
    rc = -1;
  }

  rnet.print();
 #endif // !WCS_PERF_PROF

  if (num_iter > 0u) {
   #ifdef WCS_HAS_VTUNE
    __itt_resume();
    __itt_task_begin(vtune_domain_react, __itt_null, __itt_null, vtune_handle_react);
   #endif // WCS_HAS_VTUNE
   #ifdef WCS_HAS_HPCTOOLKIT
    hpctoolkit_sampling_start();
   #endif // WCS_HAS_HPCTOOLKIT

    double t = rnet.compute_all_reaction_rates(num_iter);

   #ifdef WCS_HAS_VTUNE
    __itt_task_end(vtune_domain_react);
    __itt_pause();
   #endif // WCS_HAS_VTUNE
   #ifdef WCS_HAS_HPCTOOLKIT
    hpctoolkit_sampling_stop();
   #endif // WCS_HAS_HPCTOOLKIT

    std::cout << "Time to compute reactions rates "
              << num_iter << " times: " << t << " sec" << std::endl;
  }

  return rc;
}
