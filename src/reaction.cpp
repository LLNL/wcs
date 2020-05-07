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


#define OPTIONS "ho:"
static const struct option longopts[] = {
    {"help",    no_argument,  0, 'h'},
    {"outfile", required_argument,  0, 'o'},
    { 0, 0, 0, 0 },
};

void print_usage(const std::string exec, int code)
{
  std::cerr <<
    "Usage: " << exec << " <filename>.graphml\n"
    "    Load a reaction network graph (<filename>.graphml)\n"
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
    "    -o, --outfile\n"
    "            Specify the output file name\n"
    "\n";
  exit(code);
}

using graph_t  = boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::bidirectionalS,
    wcs::Vertex,
    wcs::Edge,
    boost::no_property,
    boost::vecS>;


void traverse(const wcs::Network& rnet)
{
  const wcs::Network::graph_t& g = rnet.graph();
  using r_prop_t = wcs::Reaction<wcs::Network::v_desc_t>;
  using s_prop_t = wcs::Species;

  std::cout << "Species:";
  for(const auto& vd : rnet.species_list()) {
    const auto& sv = g[vd];
    const auto& sp = sv.property<s_prop_t>();
    std::cout << ' ' << sv.get_label() << '[' << sp.get_count() << ']';
  }

  std::cout << "\n\nReactions:\n";
  for(const auto& vd : rnet.reaction_list()) {
    using directed_category
      = typename boost::graph_traits<wcs::Network::graph_t>::directed_category;
    constexpr bool is_bidirectional 
      = std::is_same<directed_category, boost::bidirectional_tag>::value;

    const auto& rv = g[vd];
    const auto& rp = rv.property<r_prop_t>();
    std::cout << "  " << rv.get_label()
              << " with rate constant " << rp.get_rate_constant() << " :";

    std::cout << " produces";
    for(const auto vi_out : boost::make_iterator_range(boost::out_edges(vd, g))) {
      std::cout << ' ' << g[boost::target(vi_out, g)].get_label();
    }

    if constexpr (is_bidirectional) {
      std::cout << " from a set of reactants";
      for(const auto vi_in : boost::make_iterator_range(boost::in_edges(vd, g))) {
        const auto& sv = g[boost::source(vi_in, g)];
        const auto& sp = sv.property<s_prop_t>();
        std::cout << ' ' << sv.get_label() << " [" << sp.get_count() << "]";
      }
    }

    std::cout << std::endl << "    by the rate " << rp.get_rate()
              << " <= {" << rp.get_rate_formula() << "}" << std::endl;
  }
}

int main(int argc, char** argv)
{
  int c;
  int rc = 0;
  std::string outfile;

  while ((c = getopt_long(argc, argv, OPTIONS, longopts, NULL)) != -1) {
    switch (c) {
      case 'h': /* --help */
        print_usage(argv[0], 0);
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

  traverse(rnet);
  return rc;
}
