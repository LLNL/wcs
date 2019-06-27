#include <getopt.h>
#include <boost/filesystem.hpp>
#include "reaction_network/network.hpp"
#include "utils/write_graphviz.hpp"
#include "utils/rngen.hpp"
#include <fstream>
#include <cmath>
#include <limits>
#include <vector>
#include <set>

extern "C" {
#if HAVE_CONFIG_H
#include "config.h"
#endif
}

#define OPTIONS "ho:s:i:"
static const struct option longopts[] = {
    {"help",    no_argument,  0, 'h'},
    {"outfile", required_argument,  0, 'o'},
    {"seed", required_argument,  0, 's'},
    {"iter", required_argument,  0, 'i'},
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
    "\n"
    "    -s, --seed\n"
    "            Specify the seed for random number generator. Without this it"
    "            will use a value dependent on the current system clock.\n"
    "\n"
    "    -i, --iter\n"
    "            Specify the number of reaction iterations to run.\n"
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

using rng_t = wcs::RNGen<std::uniform_int_distribution, unsigned>;
using priority_t = std::pair<wcs::etime_t, wcs::Network::v_desc_t>;
using event_queue_t = std::vector<priority_t>;

std::string show_species_names(const wcs::Network& rnet)
{
  const wcs::Network::graph_t& g = rnet.graph();

  std::string str("Species: ");
  for(const auto& vd : rnet.species_list()) {
    const auto& sv = g[vd]; // vertex (property) of the species
    str += '\t' + sv.get_label();
  }
  return str;
}

std::string show_reaction_names(const wcs::Network& rnet)
{
  const wcs::Network::graph_t& g = rnet.graph();

  std::string str("Reaction:");
  for(const auto& vd : rnet.reaction_list()) {
    const auto& rv = g[vd]; // vertex (property) of the reaction
    str += '\t' + rv.get_label();
  }
  return str;
}


std::string show_species_counts(const wcs::Network& rnet, wcs::etime_t t)
{
  using s_prop_t = wcs::Species;
  const wcs::Network::graph_t& g = rnet.graph();

  std::string str = std::to_string(t);
  for(const auto& vd : rnet.species_list()) {
    const auto& sv = g[vd]; // vertex (property) of the species
    const auto& sp = sv.property<s_prop_t>(); // detailed vertex property data of the species
    str += '\t' + std::to_string(sp.get_count());
  }
  return str;
}

std::string show_reaction_rates(const wcs::Network& rnet, wcs::etime_t t)
{
  using r_prop_t = wcs::Reaction<wcs::Network::v_desc_t>;
  const wcs::Network::graph_t& g = rnet.graph();

  std::string str = std::to_string(t);
  for(const auto& vd : rnet.reaction_list()) {
    const auto& rv = g[vd]; // vertex (property) of the reaction
    const auto& rp = rv.property<r_prop_t>(); // detailed vertex property data of the reaction
    str += '\t' + std::to_string(rp.get_rate());
  }
  return str;
}

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
std::pair<priority_t, bool> fire_reaction(const wcs::Network& rnet, event_queue_t& heap, rng_t& rgen)
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
  const wcs::etime_t t_firing = firing.first; // time to fire the reaction
  if ((t_firing == std::numeric_limits<wcs::etime_t>::infinity()) ||
      (t_firing >= wcs::Network::get_etime_ulimit())) {
    return std::make_pair(firing, false);
  } 

  // species to update as a result of the reaction fired
  std::vector<v_desc_t> updating_species;
  std::set<v_desc_t> affected_reactions;

  // product species
  for(const auto ei_out : boost::make_iterator_range(boost::out_edges(vd_firing, g))) {
    const auto vd_updating = boost::target(ei_out, g);
    updating_species.push_back(vd_updating);
    const auto& sv_updating = g[vd_updating];
    auto& sp_updating = sv_updating.property<s_prop_t>();
    const auto stoichio = g[ei_out].get_stoichiometry_ratio();
    if (!sp_updating.inc_count(stoichio)) {
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
    auto& sp_updating = sv_updating.property<s_prop_t>();
    const auto stoichio = g[ei_in].get_stoichiometry_ratio();
    if (!sp_updating.dec_count(stoichio)) {
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
    next_time = (rate_new <= static_cast<wcs::etime_t>(0))?
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
    if (rate_new <= static_cast<wcs::etime_t>(0)) {
      new_time = wcs::Network::get_etime_ulimit();
    } else if (rate_old <= static_cast<wcs::etime_t>(0)) {
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

int main(int argc, char** argv)
{
  int c;
  int rc = 0;
  unsigned seed = 0;
  unsigned num_iter = 10;
  std::string outfile;

  while ((c = getopt_long(argc, argv, OPTIONS, longopts, NULL)) != -1) {
    switch (c) {
      case 'h': /* --help */
        print_usage(argv[0], 0);
        break;
      case 'o': /* --outfile */
        outfile = std::string(optarg);
        break;
      case 's': /* --seed */
        seed = static_cast<unsigned>(atoi(optarg));
        break;
      case 'i': /* --iter */
        num_iter = static_cast<unsigned>(atoi(optarg));
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

  event_queue_t heap;
  rng_t rgen;
  if (seed == 0u)
    rgen.set_seed();
  else
    rgen.set_seed(seed);

  constexpr unsigned uint_max = std::numeric_limits<unsigned>::max();
  rgen.param(typename rng_t::param_type(10000, uint_max-10000));

  std::cout << show_species_names(rnet) << '\t';
  std::cout << show_reaction_names(rnet) << std::endl;

  build_heap(rnet, heap, rgen);

  for (unsigned int i = 0u; i < num_iter; ++i) {
    auto p = fire_reaction(rnet, heap, rgen);
    if (!p.second) {
      break;
    }

    std::cout << show_species_counts(rnet, p.first.first) << '\t';
    std::cout << show_reaction_rates(rnet, p.first.first) << std::endl;
  }
  return rc;
}
