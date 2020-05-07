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
#include "utils/graph_factory.hpp"
#include "utils/write_graphviz.hpp"
#include <fstream>

extern "C" {
#if HAVE_CONFIG_H
#include "config.h"
#endif
}

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

  wcs::GraphFactory gfactory;
  std::string fn(argv[optind]);

  if (!gfactory.read_graphml(fn)) {
    std::cerr << "Failed to read " << fn << std::endl;
    rc = -1;
    return rc;
  }

  using graph_t  = boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::bidirectionalS,
    wcs::Vertex,
    wcs::Edge,
    boost::no_property,
    boost::vecS>;

  std::cerr << "copying graph" << std::endl;
  graph_t g_polymorphic;
  gfactory.copy_to(g_polymorphic);

  if (outfile.empty()) {
    boost::filesystem::path path = fn;
    std::string base = path.stem().string();
    outfile = base + "_poly.dot";;;
  }

  if (!wcs::write_graphviz(outfile, g_polymorphic)) {
    std::cerr << "Failed to write " << outfile << std::endl;
    rc = -1;
  }
  return rc;
}
