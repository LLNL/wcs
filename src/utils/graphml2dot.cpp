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
#include "utils/graph_factory.hpp"
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
    "    Convert a reaction network graph (<filename>.graphml)\n"
    "    to AT&T GraphViz format (<filename>.dot). The output\n"
    "    is written to a file that has the same name as the input\n"
    "    file except the extention '.dot' unless --outfile is given.\n"
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

  if (outfile.empty()) {
    boost::filesystem::path path = fn;
    std::string base = path.stem().string();
    outfile = base + ".dot";
  }

  if (!gfactory.read_graphml(fn)) {
    std::cerr << "Failed to read " << fn << std::endl;
    rc = -1;
    return rc;
  }

  if (!wcs::write_graphviz(outfile, gfactory.graph())) {
    std::cerr << "Failed to write " << outfile << std::endl;
    rc = -1;
  }

  return rc;
}
