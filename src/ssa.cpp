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

#define OPTIONS "ho:"
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

  std::string str("name:");
  for(const auto& vd : rnet.species_list()) {
    const auto& sv = g[vd];
    str += '\t' + sv.get_label();
  }
  return str;
}


std::string show_species_counts(const wcs::Network& rnet, wcs::etime_t t)
{
  const wcs::Network::graph_t& g = rnet.graph();
  using s_prop_t = wcs::Species;

  std::string str = std::to_string(t);
  for(const auto& vd : rnet.species_list()) {
    const auto& sv = g[vd];
    const auto& sp = sv.property<s_prop_t>();
    str += '\t' + std::to_string(sp.get_count());
  }
  return str;
}

/**
 * Defines the priority queue ordering by the event time of entries
 * (in ascending order).
 */
bool sooner(const priority_t& v1, const priority_t& v2) {
  return (v1.first < v2.first);
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
    const auto t = log(unsigned_max/rgen())/rate;
    heap.emplace_back(priority_t(t, vd));
  }
  std::make_heap(heap.begin(), heap.end(), sooner);
}

/**
 * Pick the reaction with the earlies time to occur, execute it, update
 * the species population, recompute the reaction rates of those affected,
 * which are linked with updating species. This follows the next reaction
 * meothod procedure.
 */
priority_t fire_reaction(const wcs::Network& rnet, event_queue_t& heap, rng_t& rgen)
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
    return priority_t();
  }

  const wcs::Network::graph_t& g = rnet.graph();

  const auto firing = heap.front();
  const wcs::etime_t t_firing = firing.first; // time to fire the reaction
  const v_desc_t& vd_firing = firing.second; // vertex descriptor of the firing reaction

#if 1
  const auto& rv_firing = g[vd_firing]; // vertex (property) of the firing reaction
  const auto& rp_firing = rv_firing.property<r_prop_t>(); // detailed vertex property of the firing reaction
  std::cout << " firing " << rv_firing.get_label() << " '"
            << rp_firing.get_rate_formula() << "' at " << t_firing << std::endl;
#endif

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
    heap.front().first = log(unsigned_max/rgen())/rate_new;
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
    e.first =  e.first * rate_old / rate_new;
  }

  std::make_heap(heap.begin(), heap.end(), sooner);
  return firing;
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

  constexpr unsigned unsigned_max = std::numeric_limits<unsigned>::max();
  rgen.param(typename rng_t::param_type(10000, unsigned_max-10000));

  std::cout << show_species_names(rnet) << std::endl;
  std::cout << show_species_counts(rnet, 0.0) << std::endl;;
  build_heap(rnet, heap, rgen);

  for (unsigned int i = 0u; i < num_iter; ++i) {
    auto p = fire_reaction(rnet, heap, rgen);
    // BGL vertex descriptor of the firing reaction
    const auto vd_firing = p.second;
    // vertex (property) of the firing reaction
    const auto& rv_firing = g[vd_firing];
    using r_prop_t = wcs::Reaction<wcs::Network::v_desc_t>;
    // detailed vertex property data of the firing reaction
    const auto& rp_firing = rv_firing.property<r_prop_t>();
     std::cout << " firing " << rv_firing.get_label() << " '"
               << rp_firing.get_rate_formula() << "' at " << p.first << std::endl;

    std::cout << show_species_counts(rnet, p.first) << std::endl;
  }
  return rc;
}
