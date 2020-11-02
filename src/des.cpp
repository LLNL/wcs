/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#include <getopt.h>
#include <fstream>
#include <cmath>
#include <limits>
#include <vector>
#include <string>
#include <boost/filesystem.hpp>
#include "reaction_network/network.hpp"
#include "utils/write_graphviz.hpp"
#include "utils/timer.hpp"
#include "sim_methods/ssa_nrm.hpp"
#include "sim_methods/ssa_direct.hpp"
#include "sim_methods/ssa_sod.hpp"


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

#if defined(WCS_HAS_ROSS)
#include "sim_methods/ssa_nrm.hpp"
#include "ross.h"
#include <cstdio>
#include <cstring>


struct wcs_state
{
   wcs::SSA_NRM* ssa;
   std::shared_ptr<wcs::Network> net_ptr;
};

typedef enum {
  WCS_INIT,
  TENTATIVE_REACTION,
} wcs_event_type;

typedef struct {
   wcs::Network::v_desc_t fired_reaction;
} wcs_message;

void wcs_init(wcs_state *s, tw_lp *lp);
void wcs_event(wcs_state *s, tw_bf *bf, wcs_message *msg, tw_lp *lp);
void wcs_event_reverse(wcs_state *s, tw_bf *bf, wcs_message *msg, tw_lp *lp);
void wcs_event_commit(wcs_state *s, tw_bf *bf, wcs_message *msg, tw_lp *lp);
void wcs_final(wcs_state *s, tw_lp *lp);

tw_lptype wcs_LPs[] = {
  {
    (init_f) wcs_init,
    (pre_run_f) NULL,
    (event_f) wcs_event,
    (revent_f) wcs_event_reverse,
    (commit_f) wcs_event_commit,
    (final_f) wcs_final,
    (map_f) NULL,
    sizeof(wcs_state)
  },
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0ul },
};

const tw_optdef wcs_opts[] = {
  TWOPT_GROUP("ROSS Model"),
  TWOPT_END(),
};


int main (int argc, char* argv[])
{
  int rc = EXIT_SUCCESS;

  tw_opt_add(wcs_opts);
  tw_init(&argc, &argv);

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

  int num_LPs_per_pe = 1;

  printf("tw_nnodes %d\n", tw_nnodes());
  printf("num_LPs_per_pe %d\n", num_LPs_per_pe);
  printf("message size %lu\n", sizeof(wcs_message));

  tw_define_lps(num_LPs_per_pe, sizeof(wcs_state));

  g_tw_lp_types = wcs_LPs;
  tw_lp_setup_types();
  tw_end();

  double t_start = wcs::get_time();
  std::cout << "Wall clock time to run simulation: "
            << wcs::get_time() - t_start << " (sec)" << std::endl;

  if (cfg.tracing || cfg.sampling) {
    ssa->finalize_recording();
  } else {
    std::cout << "Species   : " << rnet.show_species_labels("") << std::endl;
    std::cout << "FinalState: " << rnet.show_species_counts() << std::endl;
  }

  delete ssa;

  return rc;
}

void wcs_init(wcs_state *s, tw_lp *lp)
{
   s->ssa->init_des(s->net_ptr, lp->gid+147); //TODO: fill with seed from command line
}

void wcs_event(wcs_state *s, tw_bf *bf, wcs_message *msg, tw_lp *lp)
{
   wcs::SSA_NRM::priority_t firing;
   firing.first = tw_now(lp);
   firing.second = msg->fired_reaction;
   s->ssa->forward_des(firing);

   bool done=false; //TODO, put end conditions in here.
   if (!done)
   {
      const auto new_firing = s->ssa->choose_reaction();

      tw_event* next_event = tw_event_new(lp->gid, new_firing.first, lp);
      wcs_message* next_msg = reinterpret_cast<wcs_message*>(tw_event_data(next_event));
      next_msg->fired_reaction = new_firing.second;
      tw_event_send(next_event);
   }
}

void wcs_event_reverse(wcs_state *s, tw_bf *bf, wcs_message *msg, tw_lp *lp)
{
   wcs::SSA_NRM::priority_t firing;
   firing.first = tw_now(lp);
   firing.second = msg->fired_reaction;
   s->ssa->backward_des(firing);
}

void wcs_event_commit(wcs_state* s, tw_bf *bf, wcs_message *msg, tw_lp *lp)
{
   s->ssa->commit_des();
}

void wcs_final(wcs_state *s, tw_lp *lp) 
{
   
}
#endif // defined(WCS_HAS_ROSS)
