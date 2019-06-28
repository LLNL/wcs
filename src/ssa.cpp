#include <getopt.h>
#include <fstream>
#include <cmath>
#include <limits>
#include <vector>
#include <set>
#include <boost/filesystem.hpp>
#include "reaction_network/network.hpp"
#include "utils/write_graphviz.hpp"
#include "utils/rngen.hpp"
#include "utils/trace_ssa.hpp"
#include "sim_methods/ssa_nrm.hpp"

extern "C" {
#if HAVE_CONFIG_H
#include "config.h"
#endif
}

#define OPTIONS "g:hi:o:s:"
static const struct option longopts[] = {
    {"graphviz", required_argument,  0, 'g'},
    {"help",     no_argument,        0, 'h'},
    {"iter",     required_argument,  0, 'i'},
    {"outfile",  required_argument,  0, 'o'},
    {"seed",     required_argument,  0, 's'},
    { 0, 0, 0, 0 },
};

struct Config {
  Config() : seed(0u), num_iter(10u) {}

  void getopt(int& argc, char** &argv);
  void print_usage(const std::string exec, int code);

  unsigned seed;
  unsigned num_iter;

  std::string infile;
  std::string outfile;
  std::string gvizfile;
};

void Config::getopt(int& argc, char** &argv)
{
  int c;

  while ((c = getopt_long(argc, argv, OPTIONS, longopts, NULL)) != -1) {
    switch (c) {
      case 'g': /* --graphviz */
        gvizfile = std::string(optarg);
        break;
      case 'h': /* --help */
        print_usage(argv[0], 0);
        break;
      case 'i': /* --iter */
        num_iter = static_cast<unsigned>(atoi(optarg));
        break;
      case 'o': /* --outfile */
        outfile = std::string(optarg);
        break;
      case 's': /* --seed */
        seed = static_cast<unsigned>(atoi(optarg));
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
    "    -g, --graphviz\n"
    "            Specify the name of the file to export the reaction\n"
    "            network into in the GraphViz format.\n"
    "\n"
    "    -o, --outfile\n"
    "            Specify the output file name\n"
    "\n"
    "    -s, --seed\n"
    "            Specify the seed for random number generator. Without this,\n"
    "            it will use a value dependent on the current system clock.\n"
    "\n"
    "    -i, --iter\n"
    "            Specify the number of reaction iterations to run.\n"
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

  wcs::SSA_NRM ssa;
  ssa.init(rnet_ptr, cfg.num_iter, cfg.seed, true);
  ssa.run();

  if (!cfg.outfile.empty()) {
    ssa.trace().write(cfg.outfile);
  }

  return rc;
}
