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
#include "sim_methods/ssa_nrm.hpp"


#define OPTIONS "dg:hi:o:s:t:"
static const struct option longopts[] = {
    {"diag",     no_argument,        0, 'd'},
    {"graphviz", required_argument,  0, 'g'},
    {"help",     no_argument,        0, 'h'},
    {"iter",     required_argument,  0, 'i'},
    {"outfile",  required_argument,  0, 'o'},
    {"seed",     required_argument,  0, 's'},
    {"time",     required_argument,  0, 't'},
    { 0, 0, 0, 0 },
};

struct Config {
  Config()
  : seed(0u), max_iter(10u),
    max_time(wcs::Network::get_etime_ulimit()),
    tracing(false)
  {}

  void getopt(int& argc, char** &argv);
  void print_usage(const std::string exec, int code);

  unsigned seed;
  unsigned max_iter;
  wcs::sim_time_t max_time;
  bool tracing;

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
        break;
      case 'g': /* --graphviz */
        gvizfile = std::string(optarg);
        break;
      case 'h': /* --help */
        print_usage(argv[0], 0);
        break;
      case 'i': /* --iter */
        max_iter = static_cast<unsigned>(atoi(optarg));
        is_iter_set = true;
        break;
      case 'o': /* --outfile */
        outfile = std::string(optarg);
        tracing = true;
        break;
      case 's': /* --seed */
        seed = static_cast<unsigned>(atoi(optarg));
        break;
      case 't': /* --time */
        max_time = static_cast<wcs::sim_time_t>(std::stod(optarg));
        is_time_set = true;
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
    "            Specify the tracing output file name. This enables tracing.\n"
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
  ssa.init(rnet_ptr, cfg.max_iter, cfg.max_time, cfg.seed, cfg.tracing);
  ssa.run();

  if (!cfg.outfile.empty() && cfg.tracing) {
    ssa.trace().write(cfg.outfile);
  } else if (cfg.tracing) {
    ssa.trace().write(std::cout);
  }
  if (!cfg.tracing) {
    std::cout << "Species   : " << rnet.show_species_labels("") << std::endl;
    std::cout << "FinalState: " << rnet.show_species_counts() << std::endl;
  }

  return rc;
}
