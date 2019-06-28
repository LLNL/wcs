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

using rng_t = wcs::RNGen<std::uniform_int_distribution, unsigned>;
using priority_t = std::pair<wcs::sim_time_t, wcs::Network::v_desc_t>;
using event_queue_t = std::vector<priority_t>;


/**
 * Defines the priority queue ordering by the event time of entries
 * (in ascending order).
 */
bool later(const priority_t& v1, const priority_t& v2) {
  return (v1.first >= v2.first);
}

/**
 * Initialize the priority queue.
 * For each reaction, compute the time to reaction using the random number
 * generator rgen.
 */
void build_heap(const wcs::Network& rnet, event_queue_t& heap, rng_t& rgen)
{
  using r_prop_t = wcs::Reaction<wcs::Network::v_desc_t>;

  using directed_category = typename boost::graph_traits<wcs::Network::graph_t>::directed_category;
  constexpr bool is_bidirectional
    = std::is_same<directed_category, boost::bidirectional_tag>::value;

  if constexpr (!is_bidirectional) {
    WCS_THROW("Cannot get species population without in-edges.");
  }

  const wcs::Network::graph_t& g = rnet.graph();

  heap.clear();
  heap.reserve(rnet.get_num_reactions()+10);
  constexpr double unsigned_max = static_cast<double>(std::numeric_limits<unsigned>::max());

  // For each reaction, compute the time-to-reaction
  for(const auto& vd : rnet.reaction_list()) {
    const auto& rv = g[vd]; // reaction vertex
    const auto& rp = rv.property<r_prop_t>(); // detailed vertex property data
    const auto rate = rp.get_rate(); // reaction rate
    const auto rn = unsigned_max/rgen();
    const auto t = log(rn)/rate;
    heap.emplace_back(priority_t(t, vd));
  }
  std::make_heap(heap.begin(), heap.end(), later);
}

/**
 * Pick the reaction with the earlies time to occur, execute it, update
 * the species population, recompute the reaction rates of those affected,
 * which are linked with updating species. This follows the next reaction
 * meothod procedure.
 */
std::pair<priority_t, bool>
fire_reaction(const wcs::Network& rnet,
              event_queue_t& heap,
              rng_t& rgen,
              wcs::Trace& trace)
{
  using v_desc_t = wcs::Network::v_desc_t;
  using r_prop_t = wcs::Reaction<wcs::Network::v_desc_t>;
  using s_prop_t = wcs::Species;
  constexpr double unsigned_max = static_cast<double>(std::numeric_limits<unsigned>::max());


  using directed_category = typename boost::graph_traits<wcs::Network::graph_t>::directed_category;
  constexpr bool is_bidirectional
    = std::is_same<directed_category, boost::bidirectional_tag>::value;

  if constexpr (!is_bidirectional) {
    WCS_THROW("Cannot get species population without in-edges.");
  }

  if (heap.empty()) {
    return std::make_pair(priority_t(), false);
  }

  const wcs::Network::graph_t& g = rnet.graph();

  const auto firing = heap.front();
  const v_desc_t& vd_firing = firing.second; // vertex descriptor of the firing reaction
  const wcs::sim_time_t t_firing = firing.first; // time to fire the reaction
  if ((t_firing == std::numeric_limits<wcs::sim_time_t>::infinity()) ||
      (t_firing >= wcs::Network::get_etime_ulimit())) {
    return std::make_pair(firing, false);
  } 
  trace.record_reaction(t_firing, vd_firing);

  // species to update as a result of the reaction fired
  std::vector<v_desc_t> updating_species;
  std::set<v_desc_t> affected_reactions;

  // product species
  for(const auto ei_out : boost::make_iterator_range(boost::out_edges(vd_firing, g))) {
    const auto vd_updating = boost::target(ei_out, g);
    updating_species.push_back(vd_updating);
    const auto& sv_updating = g[vd_updating];
    if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
      // in case that there are other type of vertices than species or reaction
      if (sv_updating.get_type() != wcs::Vertex::_species_) continue;
    }

    auto& sp_updating = sv_updating.property<s_prop_t>();
    const auto stoichio = g[ei_out].get_stoichiometry_ratio();
    if (!sp_updating.inc_count(stoichio)) { // State update
      // TODO: For reversible implementation, this has to be carefully handled.
      // e.g., record which one was successful or not successful.
      continue;
    }

    for(const auto vi_affected : boost::make_iterator_range(boost::out_edges(vd_updating, g))) {
      const auto vd_affected = boost::target(vi_affected, g);
      if (vd_affected == vd_firing) continue;
      affected_reactions.insert(vd_affected);
    }
  }

  // reactant species
  for(const auto ei_in : boost::make_iterator_range(boost::in_edges(vd_firing, g))) {
    const auto vd_updating = boost::source(ei_in, g);
    // TODO: in case that there are other type of vertices than species or reaction
    // need to check if the vertex is species type
    updating_species.push_back(vd_updating);
    const auto& sv_updating = g[vd_updating];
    if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
      // in case that there are other type of vertices than species or reaction
      if (sv_updating.get_type() != wcs::Vertex::_species_) continue;
    }

    auto& sp_updating = sv_updating.property<s_prop_t>();
    const auto stoichio = g[ei_in].get_stoichiometry_ratio();
    if (!sp_updating.dec_count(stoichio)) { // State update
      // TODO: For reversible implementation, this has to be carefully handled.
      // e.g., record which one was successful or not successful.
      continue;
    }

    for(const auto vi_affected : boost::make_iterator_range(boost::out_edges(vd_updating, g))) {
      const auto vd_affected = boost::target(vi_affected, g);
      if (vd_affected == vd_firing) continue;
      affected_reactions.insert(vd_affected);
    }
  }

  { // process the firing reaction at the end of the heap
    const auto rate_new = rnet.set_reaction_rate(vd_firing);
    const auto rn = unsigned_max/rgen();
    auto& next_time = heap.front().first;
    next_time = (rate_new <= static_cast<wcs::sim_time_t>(0))?
                   wcs::Network::get_etime_ulimit() :
                   log(rn)/rate_new;
    if (isnan(next_time)) {
      next_time = wcs::Network::get_etime_ulimit();
    }
  }

  // update the event time of the rest of affected reactions
  for(auto& e: heap) {
    // TODO: need a better facility to locate the affected reaction entry in
    // the heap. red-black tree perhaps.
    if (affected_reactions.count(e.second) == 0u) {
      continue;
    }
    const auto& rv_affected = g[e.second];
    auto& rp_affected = rv_affected.property<r_prop_t>();
    const auto rate_old = rp_affected.get_rate();
    const auto rate_new = rnet.set_reaction_rate(e.second);
    auto& new_time = e.first;
    if (rate_new <= static_cast<wcs::sim_time_t>(0)) {
      new_time = wcs::Network::get_etime_ulimit();
    } else if (rate_old <= static_cast<wcs::sim_time_t>(0)) {
      const auto rn = unsigned_max/rgen();
      new_time = log(rn)/rate_new;
    } else {
      new_time = new_time * rate_old / rate_new;
    }
    if (isnan(new_time)) {
      new_time = wcs::Network::get_etime_ulimit();
    }
  }

  std::make_heap(heap.begin(), heap.end(), later);
  return std::make_pair(firing, true);
}


void run_SSA_next_reaction_method(
        wcs::Network& rnet,
        wcs::Trace& trace,
        rng_t& rgen,
        unsigned int max_iter)
{
  wcs::sim_time_t sim_time = 0.0;
  std::cout << rnet.show_species_labels() << '\t';
  std::cout << rnet.show_reaction_labels() << std::endl;
  std::cout << sim_time << rnet.show_species_counts() << '\t';
  std::cout << sim_time << rnet.show_reaction_rates() << std::endl;

  event_queue_t heap;
  build_heap(rnet, heap, rgen);

  for (unsigned int i = 0u; i < max_iter; ++i) {
    auto p = fire_reaction(rnet, heap, rgen, trace);
    if (!p.second) {
      break;
    }
    sim_time += p.first.first;
    std::cout << sim_time << rnet.show_species_counts() << '\t';
    std::cout << sim_time << rnet.show_reaction_rates() << std::endl;
  }
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

  wcs::TraceSSA trace;
  trace.record_initial_condition(rnet_ptr);

  if (!cfg.gvizfile.empty() &&
      !wcs::write_graphviz(cfg.gvizfile, g))
  {
    std::cerr << "Failed to write " << cfg.gvizfile << std::endl;
    rc = -1;
  }

  rng_t rgen;
  if (cfg.seed == 0u)
    rgen.set_seed();
  else
    rgen.set_seed(cfg.seed);

  constexpr unsigned uint_max = std::numeric_limits<unsigned>::max();
  rgen.param(typename rng_t::param_type(1000, uint_max-1000));

  run_SSA_next_reaction_method(rnet, trace, rgen, cfg.num_iter);

  if (!cfg.outfile.empty()) {
    trace.write(cfg.outfile);
  }

  return rc;
}
