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

#include "reaction_network/network.hpp"
#include "utils/graph_factory.hpp"
#include "utils/input_filetype.hpp"
#include "utils/generate_cxx_code.hpp"
#include "utils/timer.hpp"
#include <type_traits> // is_same<>
#include <algorithm> // lexicographical_compare(), sort()
#include <limits> // numeric_limits

#if defined(WCS_HAS_SBML)
#include <sbml/SBMLTypes.h>
#include <sbml/common/extern.h>
#endif // defined(WCS_HAS_SBML)

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  @{ */

sim_time_t Network::m_etime_ulimit = std::numeric_limits<sim_time_t>::infinity();

void Network::load(const std::string filename, const bool reuse)
{
  input_filetype fn(filename);
  input_filetype::input_type filetype = fn.detect();
  if (filetype == input_filetype::input_type::_graphml_) {
    loadGraphML(filename);
  } else if (filetype == input_filetype::input_type::_sbml_) {
    loadSBML(filename, reuse);
  } else if (filetype == input_filetype::input_type::_ioerror_) {
    WCS_THROW("Could not find the requested file.");
    return;
  } else if (filetype == input_filetype::input_type::_unknown_) {
    WCS_THROW("Unknown filetype. Please select a reaction network graph "
              "(<filename>.graphml) or a Systems Biology Markup Language "
              "(SBML) file (<filename>.xml)\n");
    return;
  }
}

void Network::loadGraphML(const std::string graphml_filename)
{
  #if !defined(WCS_HAS_EXPRTK)
  WCS_THROW("Must enable ExprTk for .graphml files.");
  return;
  #endif
  ::wcs::GraphFactory gfactory;

  if (!gfactory.read_graphml(graphml_filename)) {
    WCS_THROW("Failed to read " + graphml_filename);
    return;
  }
  gfactory.copy_to(m_graph);
}

void Network::print_parameters_of_reactions(
  const params_map_t& dep_params_f,
  const params_map_t& dep_params_nf,
  const rate_rules_dep_t& rate_rules_dep_map)
{
  #if !defined(WCS_HAS_EXPRTK)
  typename params_map_t::const_iterator pit;
  typename rate_rules_dep_t::const_iterator rrdit;
  using std::operator<<;
  for(auto &e : dep_params_f ) {
    std::string reaction_name = e.first;
    std::vector<std::string> params_f = e.second;
    std::cout << "Reaction: " << reaction_name << " Params: ";
    for (auto& x: params_f) {
      std::cout << x << ", ";
      rrdit = rate_rules_dep_map.find(x);
      if (rrdit != rate_rules_dep_map.cend()) {
        std::cout << " (";
        std::set<std::string> params_rr = rrdit->second;
        for (auto& w: params_rr) {
          std::cout << w << ", ";
        }
        std::cout << ") ";
      }
    }
    std::cout << "/ ";
    pit = dep_params_nf.find(reaction_name);
    if (pit != dep_params_nf.cend()) {
      std::vector<std::string> params_nf = pit->second;
      for (auto& x: params_nf) {
        std::cout << x << ", ";
        rrdit = rate_rules_dep_map.find(x);
        if (rrdit != rate_rules_dep_map.cend()) {
          std::cout << " (";
          std::set<std::string> params_rr = rrdit->second;
          for (auto& w: params_rr) {
            std::cout << w << ", ";
          }
          std::cout << ") ";
        }
      }
    }
    std::cout << std::endl;
  }
  #endif // !defined(WCS_HAS_EXPRTK
}

void Network::loadSBML(const std::string sbml_filename, const bool reuse)
{
  #if defined(WCS_HAS_SBML)

  ::wcs::GraphFactory gfactory;
  LIBSBML_CPP_NAMESPACE::SBMLReader reader;
  const LIBSBML_CPP_NAMESPACE::SBMLDocument* document
    = reader.readSBML(sbml_filename);

  if (document == nullptr) {
    WCS_THROW("Failed to read the SBML file " + sbml_filename);
    return;
  }

  const LIBSBML_CPP_NAMESPACE::Model* model = document->getModel();

  const unsigned num_errors = document->getNumErrors();

  if (num_errors > 0u) {
    document->printErrors(std::cerr);

    delete document;
    WCS_THROW(std::to_string(num_errors) + " error(s) in reading " +
              sbml_filename);
    return;
  }

  if (model == nullptr) {
    delete document;
    WCS_THROW("Failed to get model from " + sbml_filename);
    return;
  }

  #if !defined(WCS_HAS_EXPRTK)

  const std::string lib_filename = get_libname_from_model(sbml_filename);
  generate_cxx_code code_generator(lib_filename, !reuse);

  typename params_map_t::const_iterator pit;
  struct timespec t_1, t_2; 
  code_generator.generate_code(*model,
        m_dep_params_f, m_dep_params_nf, m_rate_rules_dep_map);
  const std::string library_file = code_generator.compile_code();

  using std::operator<<;
  std::cerr << "Constructing a graph from the SBML model ..." << std::endl;
  // print_parameters_of_reactions(m_dep_params_f, m_dep_params_nf,
  //                               m_rate_rules_dep_map);

  gfactory.convert_to(*model, m_graph, library_file,
                      m_dep_params_f, m_dep_params_nf,
                      m_rate_rules_dep_map);

  #else
  gfactory.convert_to(*model, m_graph, "",{},{},{});
  #endif // !defined(WCS_HAS_EXPRTK)

  delete document;

  #else
  WCS_THROW("SBML is not enabled.");
  #endif // defined(WCS_HAS_SBML)
  std::cerr << "The SBML model has been loaded into a graph!" << std::endl;

}

void Network::init()
{
  #if !defined(WCS_HAS_EXPRTK) && !defined(WCS_HAS_SBML)
  WCS_THROW("Must enable either ExprTk or SBML.");
  #endif
  const size_t num_vertices = get_num_vertices();

  m_reactions.reserve(num_vertices);
  m_species.reserve(num_vertices);

  v_iter_t vi, vi_end;

  for (boost::tie(vi, vi_end) = boost::vertices(m_graph); vi != vi_end; ++vi) {
    const v_prop_t& v = m_graph[*vi];
    const auto vt = static_cast<v_prop_t::vertex_type>(v.get_typeid());
    if (vt == v_prop_t::_species_) {
      m_species.emplace_back(*vi);
    } else {
      using directed_category = boost::graph_traits<graph_t>::directed_category;
      constexpr bool is_bidirectional
        = std::is_same<directed_category, boost::bidirectional_tag>::value;

      const v_desc_t reaction = *vi;
      // `involved_species` include all the species involved in the reaction:
      // the rate-determining species as the reactants, the enzymes and the
      // inhibiters, as well as the products.
      // They may exclude products depending on how formula parsing is
      // implemented.
      s_involved_t involved_species;

      s_involved_t products;

      #if !defined(WCS_HAS_EXPRTK)
      typename params_map_t::const_iterator pit, pit_nf;
      std::vector<std::string> params_reactants;
      std::string reaction_name = m_graph[*vi].get_label();
      pit = m_dep_params_f.find(reaction_name);
      if (pit != m_dep_params_f.end()) {
        params_reactants=pit->second;
      }
      #endif // !defined(WCS_HAS_EXPRTK)

      if constexpr (is_bidirectional) {
        for(const auto ei_in :
            boost::make_iterator_range(boost::in_edges(reaction, m_graph))) {
          v_desc_t reactant = boost::source(ei_in, m_graph);

          #if !defined(WCS_HAS_EXPRTK)
          //check for the reactants which are not actually reactants and put their stoichiometry 0
          if ( std::find(params_reactants.begin(), params_reactants.end(),
          m_graph[reactant].get_label()) == params_reactants.end()) {
            //std::cout << "Not reactant " << m_graph[reactant].get_label() << std::endl;
            m_graph[ei_in].set_stoichiometry_ratio(0);
          }
          #endif // !defined(WCS_HAS_EXPRTK)
          const auto st = m_graph[ei_in].get_stoichiometry_ratio();
          involved_species.insert(std::make_pair(m_graph[reactant].get_label(),
                                                 std::make_pair(reactant, st)));
                                              
        }
      }

      for(const auto ei_out :
          boost::make_iterator_range(boost::out_edges(reaction, m_graph))) {
        v_desc_t product = boost::target(ei_out, m_graph);
        products.insert(std::make_pair(m_graph[product].get_label(),
                                       std::make_pair(product, 1)));
      }

      m_reactions.emplace_back(reaction);

      auto& r = m_graph[*vi].checked_property< Reaction<v_desc_t> >();

      #if !defined(WCS_HAS_EXPRTK)
      pit = m_dep_params_f.find(reaction_name);
      pit_nf = m_dep_params_nf.find(reaction_name);
      typename rate_rules_dep_t::const_iterator rrdit;

      if (pit != m_dep_params_f.cend()) {
        if (pit_nf != m_dep_params_nf.cend()) {
          const std::vector<std::string>& params_fv = pit->second;
          const std::vector<std::string>& params_nfv = pit_nf->second;
          std::vector<std::string> fparams, nfparams;
          for (const auto& x: params_fv) {
            rrdit = m_rate_rules_dep_map.find(x);
            if (rrdit != m_rate_rules_dep_map.cend()) {
              const std::set<std::string>& params_rr = rrdit->second;
              for (const auto& w: params_rr) {
                fparams.push_back(w);
              }
            } else {
              fparams.push_back(x);
            }
          }

          for (const auto& x: params_nfv) {
            rrdit = m_rate_rules_dep_map.find(x);
            if (rrdit != m_rate_rules_dep_map.cend()) {
              const std::set<std::string>& params_rr = rrdit->second;
              for (const auto& w: params_rr) {
                nfparams.push_back(w);
              }
            } else{
              nfparams.push_back(x);
            }
          }

          r.set_rate_inputs(involved_species, fparams, nfparams);
        }
      } else {
        WCS_THROW("No function with the name " + reaction_name);
      }
      #else
      r.set_rate_inputs(involved_species);
      #endif // !defined(WCS_HAS_EXPRTK)

      r.set_products(products);

      set_reaction_rate(*vi);
    }
  }

  sort_species();
  build_index_maps();

  m_pid = unassigned_partition;
}

/// Overwrite the reaction rate to a given value
void Network::set_reaction_rate(const Network::v_desc_t r,
                                const wcs::reaction_rate_t rate) const
{
  const auto& rv = m_graph[r]; // vertex (property) of the reaction
  auto& rp = rv.property<r_prop_t>(); // detailed vertex property data
  rp.set_rate(rate);
}

/**
 * Computes the reaction rate based on the population of the reaction driving
 * species and the reaction constant.
 */
reaction_rate_t Network::set_reaction_rate(const Network::v_desc_t r) const
{
  auto & rprop = m_graph[r].checked_property< Reaction<v_desc_t> >();
  const auto& ri = rprop.get_rate_inputs();
  std::vector<reaction_rate_t> params;
  // GG: rate constant is part of the Reaction object
  params.reserve(ri.size()+1u); // reserve space for species count and rate constant

  for (auto driver : ri) { // add species counts here and the rate constant will be appended later
    const auto& s = m_graph[driver.first].checked_property<Species>();
    const stoic_t num_same = driver.second;
    // A reaction may take a same reactant species multiple times.
    // e.g., X + X -> Y
    // In such a case, num_same_reactants is more than 1. Two in the example.
    // In the example, the reaction rate is computed as [X]([X]-1)/2
    // TODO: move this logic into ReactionBase::set_rate_fn()
    species_cnt_t n = s.get_count();

    if (num_same == static_cast<stoic_t>(0)) {
        // This is a parameter that the reaction rate is dependent on but
        // not modified by the reaction (i.e., neither reactant nor product)
        params.push_back(static_cast<reaction_rate_t>(n));
    } else {
      params.push_back(static_cast<reaction_rate_t>(n));
    }
  }
  return rprop.calc_rate(std::move(params));
}

double Network::compute_all_reaction_rates(const unsigned n) const
{
  double t_start = get_time();
  for (unsigned i = 0u; i < n; i++) {
    for (const auto& r: reaction_list()) {
      auto & rprop = m_graph[r].checked_property< Reaction<v_desc_t> >();
      const auto& ri = rprop.get_rate_inputs();
      std::vector<reaction_rate_t> params;
      // GG: rate constant is part of the Reaction object
      params.reserve(ri.size()+1u); // reserve space for species count and rate constant

      for (auto driver : ri) { // add species counts here and the rate constant will be appended later
        const auto& s = m_graph[driver.first].checked_property<Species>();
        //const stoic_t num_same = driver.second;
        species_cnt_t n = s.get_count();
        params.push_back(static_cast<reaction_rate_t>(n));
      }

      rprop.calc_rate(std::move(params));
    }
  }
  return get_time() - t_start;
}

reaction_rate_t Network::get_reaction_rate(const Network::v_desc_t r) const
{
  const auto& rv = m_graph[r]; // vertex (property) of the reaction
  const auto& rp = rv.property<r_prop_t>(); // detailed vertex property data
  return rp.get_rate();
}

void Network::sort_species()
{
  const auto& g = m_graph;
  std::sort(m_species.begin(), m_species.end(),
            [&g](const Network::v_desc_t lhs, const Network::v_desc_t rhs) {
              const auto lstr = g[lhs].get_label();
              const auto rstr = g[rhs].get_label();
              return std::lexicographical_compare(lstr.begin(), lstr.end(),
                                                  rstr.begin(), rstr.end());
            });
}

Network::v_desc_t Network::find_species(const std::string& label)
{
  const auto& g = m_graph;
  species_list_t::const_iterator it
    = std::lower_bound(m_species.begin(), m_species.end(), label,
            [&g](const Network::v_desc_t lhs, const std::string& l) {
              const auto n = g[lhs].get_label();
              return std::lexicographical_compare(n.begin(), n.end(),
                                                  l.begin(), l.end());
            });
  return *it;
}


size_t Network::get_num_vertices() const
{
  return m_graph.m_vertices.size();
}

size_t Network::get_num_species() const
{
  return m_species.size();
}

size_t Network::get_num_reactions() const
{
  return (get_num_vertices() - get_num_species());
}

size_t Network::get_num_vertices_of_type(const Network::vertex_type vt) const
{
  switch(vt) {
    case v_prop_t::_undefined_: return get_num_vertices();
    case v_prop_t::_species_: return get_num_species();
    case v_prop_t::_reaction_: return get_num_reactions();
    case v_prop_t::_num_vertex_types_: return 0ul;
  }
  return 0ul;
}

const Network::graph_t& Network::graph() const
{
  return m_graph;
}

const Network::reaction_list_t& Network::reaction_list() const
{
  return m_reactions;
}

const Network::species_list_t& Network::species_list() const
{
  return m_species;
}

void Network::set_etime_ulimit(const sim_time_t t)
{
  m_etime_ulimit = t;
}

sim_time_t Network::get_etime_ulimit()
{
  return m_etime_ulimit;
}

/**
 * Check if the condition for reaction is satisfied such as whether a
 * sufficient number of reactants exist. If not, the reaction rate is
 * set to zero.
 */
bool Network::check_reaction(const wcs::Network::v_desc_t r) const
{
  if (m_graph[r].get_type() != wcs::Vertex::_reaction_) {
    WCS_THROW(m_graph[r].get_label() + " is not a reaction.");
    return false;
  }

  // reactant species
  for (const auto ei_in : boost::make_iterator_range(boost::in_edges(r, m_graph))) {
    const auto vd_reactant = boost::source(ei_in, m_graph);
    const auto& sv_reactant = m_graph[vd_reactant];
    if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
      // in case that there are other type of vertices than species or reaction
      if (sv_reactant.get_type() != wcs::Vertex::_species_) continue;
    }

    const auto& sp_reactant = sv_reactant.property<Species>();
    const auto stoichio = m_graph[ei_in].get_stoichiometry_ratio();
    if (!sp_reactant.dec_check(stoichio)) {
     #if 0  // change to 1 to take all messages (0 default)
      using std::operator>>;
      std::cerr << "reaction " << m_graph[r].get_label()
                << " has insufficient amount of reactants "
                << sv_reactant.get_label() << " (" << stoichio << " < "
                << sp_reactant.get_count() << ")" << std::endl;
     #endif
      // check if reaction is possible, i.e., decrement is possible
      set_reaction_rate(r, 0.0);
      return false;
    }
  }

  // product species
  for (const auto ei_out : boost::make_iterator_range(boost::out_edges(r, m_graph))) {
    const auto vd_product = boost::target(ei_out, m_graph);
    const auto& sv_product = m_graph[vd_product];
    if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
      // in case that there are other type of vertices than species or reaction
      if (sv_product.get_type() != wcs::Vertex::_species_) continue;
    }

    const auto& sp_product = sv_product.property<Species>();
    const auto stoichio = m_graph[ei_out].get_stoichiometry_ratio();
    if (!sp_product.inc_check(stoichio)) {
      // check if reaction is possible, i.e., increment is possible
      using std::operator>>;
      std::cerr << "reaction " << m_graph[r].get_label()
                << " cannot increment the amount of products."
                << " To enable 64-bit counter, rebuild using the cmake"
                << " option '-DWCS_64BIT_CNT=ON'."
                << std::endl;
      set_reaction_rate(r, 0.0);
      return false;
    }
  }

  // Alternatively, the rate computation may involve dividing the species count
  // by the stoichiometry ratio. If the species count is less than the ratio,
  // the resultant rate would be zero.

  return true;
}

std::tuple<reaction_rate_t, reaction_rate_t, reaction_rate_t>
Network::find_min_max_rate() const
{
  reaction_rate_t r_min = std::numeric_limits<reaction_rate_t>::max();
  reaction_rate_t r_max = std::numeric_limits<reaction_rate_t>::min();
  reaction_rate_t r_sum = static_cast<reaction_rate_t>(0);

  for(const auto& vd : reaction_list()) {
    const auto& rv = m_graph[vd]; // vertex (property) of the reaction
    const auto& rp = rv.property<r_prop_t>(); // detailed vertex property data
    const auto r = rp.get_rate();
    r_min = std::min(r_min, r);
    r_max = std::max(r_max, r);
    r_sum += r;
  }

  return std::make_tuple(r_min, r_max, r_sum);
}

std::string Network::show_species_labels(const std::string title) const
{
  std::string str(title);
  str.reserve(str.size() + get_num_species()*30);

  for(const auto& vd : species_list()) {
    const auto& sv = m_graph[vd]; // vertex (property) of the species
    str += '\t' + sv.get_label();
  }
  return str;
}

std::string Network::show_reaction_labels(const std::string title) const
{
  std::string str(title);
  str.reserve(str.size() + get_num_reactions()*20);

  for(const auto& vd : reaction_list()) {
    const auto& rv = m_graph[vd]; // vertex (property) of the reaction
    str += '\t' + rv.get_label();
  }
  return str;
}

std::string Network::show_species_counts() const
{
  std::string str;
  str.reserve(get_num_species()*10);

  for(const auto& vd : species_list()) {
    const auto& sv = m_graph[vd]; // vertex (property) of the species
    const auto& sp = sv.property<Species>(); // detailed vertex property data
    str += '\t' + std::to_string(sp.get_count());
  }
  return str;
}

std::string Network::show_reaction_rates() const
{
  std::string str;
  str.reserve(get_num_reactions()*15);

  for(const auto& vd : reaction_list()) {
    const auto& rv = m_graph[vd]; // vertex (property) of the reaction
    const auto& rp = rv.property<r_prop_t>(); // detailed vertex property data
    str += '\t' + std::to_string(rp.get_rate());
  }
  return str;
}

/**
 * Build the maps from the vertex descriptor to the index for species and
 * reactions respectively.
 */
void Network::build_index_maps()
{
  if constexpr (!is_vertex_list_ordered::value) {
    std::string errmsg = std::string("Currently, the mapping from a vertex ")
                       + "index to a vertex descriptor is based on the "
                       + "assumption that the vertex index is a integral "
                       + "sequence number starting from 0 to the number of "
                       + "vertices in the graph container minus one. Aslo, "
                       + "it is assumed that the container is ordered such "
                       + "that iterating over the container results in the "
                       + "same order of vertices.";
    WCS_THROW(errmsg);
    return;
  }

  { // build a map from the vertex descriptor and the index for species
    m_s_idx_map.clear();
    m_s_idx_map.reserve(m_species.size());

    v_idx_t sidx = static_cast<v_idx_t>(0u);
    for (const auto& sd : m_species) {
      m_s_idx_map[sd] = sidx++;
    }
  }
  { // build a map from the vertex descriptor and the index for reactions
    m_r_idx_map.clear();
    m_r_idx_map.reserve(m_reactions.size());

    v_idx_t ridx = static_cast<v_idx_t>(0u);
    for (const auto& rd : m_reactions) {
      m_r_idx_map[rd] = ridx++;
    }
  }
}

const Network::map_desc2idx_t& Network::get_reaction_map() const
{
  return m_r_idx_map;
}

const Network::map_desc2idx_t& Network::get_species_map() const
{
  return m_s_idx_map;
}

v_idx_t Network::reaction_d2i(v_desc_t d) const
{
  return m_r_idx_map.at(d);
}

Network::v_desc_t Network::reaction_i2d(v_idx_t i) const
{
  return m_reactions.at(i);
}

v_idx_t Network::species_d2i(v_desc_t d) const
{
  return m_s_idx_map.at(d);
}

Network::v_desc_t Network::species_i2d(v_idx_t i) const
{
  return m_species.at(i);
}

void Network::set_partition(const map_idx2desc_t& idx2vd,
                            const std::vector<partition_id_t>& parts,
                            const partition_id_t my_pid)
{
  if (idx2vd.size() != parts.size()) {
    std::string errmsg =
      "Inconsistent sizes between the vertex map and the number of partitions!";
    WCS_THROW(errmsg);
    return;
  }
  m_pid = my_pid;

  m_my_reactions.clear();
  m_my_reactions.reserve(m_reactions.size());
  m_my_species.clear();
  m_my_species.reserve(m_species.size());

  size_t i = 0u;

  for (const auto vd : idx2vd) {
    auto& v = m_graph[vd]; // vertex (property) of the reaction
    const partition_id_t pid = parts[i++];
    v.set_partition(pid);

    if (pid == my_pid) {
      const auto vt = static_cast<v_prop_t::vertex_type>(v.get_typeid());

      if (vt == v_prop_t::_species_) {
        m_my_species.emplace_back(vd);
      } else if (vt == v_prop_t::_reaction_) {
        m_my_reactions.emplace_back(vd);
      }
    }
  }
}

void Network::set_partition(const std::vector<partition_id_t>& parts,
                            const partition_id_t my_pid)
{
  if (get_num_vertices() != parts.size()) {
    std::string errmsg =
      "Inconsistent sizes between the number of verticesand the number of partitions!";
    WCS_THROW(errmsg);
    return;
  }
  m_pid = my_pid;

  m_my_reactions.clear();
  m_my_reactions.reserve(m_reactions.size());
  m_my_species.clear();
  m_my_species.reserve(m_species.size());

  size_t i = 0u;

  v_iter_t vi, vi_end;

  for (boost::tie(vi, vi_end) = boost::vertices(m_graph); vi != vi_end; ++vi) {
    auto& v = m_graph[*vi]; // vertex (property) of the reaction
    const partition_id_t pid = parts[i++];
    v.set_partition(pid);

    if (pid == my_pid) {
      const auto vt = static_cast<v_prop_t::vertex_type>(v.get_typeid());

      if (vt == v_prop_t::_species_) {
        m_my_species.emplace_back(*vi);
      } else if (vt == v_prop_t::_reaction_) {
        m_my_reactions.emplace_back(*vi);
      }
    }
    // TODO: else if it is not connected to any local vertex
    // deallocate the proporty specific to the vertex type
  }
}

const Network::reaction_list_t& Network::my_reaction_list() const
{
  return m_my_reactions;
}

const Network::reaction_list_t& Network::my_species_list() const
{
  return m_my_species;
}

partition_id_t Network::get_partition_id() const
{
  return m_pid;
}

void Network::print() const
{
  using s_prop_t = wcs::Species;

  std::cout << "Species:";
  for(const auto& vd : species_list()) {
    const auto& sv = m_graph[vd];
    const auto& sp = sv.property<s_prop_t>();
    std::cout << ' ' << sv.get_label() << '[' << sp.get_count() << ']';
  }
  size_t num_inactive = 0ul;

  std::cout << "\n\nReactions:\n";
  for(const auto& vd : reaction_list()) {
    using directed_category
      = typename boost::graph_traits<wcs::Network::graph_t>::directed_category;
    constexpr bool is_bidirectional
      = std::is_same<directed_category, boost::bidirectional_tag>::value;

    const auto& rv = m_graph[vd]; // reaction vertex
    const auto& rp = rv.property<r_prop_t>(); // reaction vertex property
    std::cout << "  " << rv.get_label()
              << " with rate constant " << rp.get_rate_constant() << " :";

    bool inactive = false;

    std::cout << " produces";
    for(const auto vi_out : boost::make_iterator_range(boost::out_edges(vd, m_graph))) {
      const auto& sv = m_graph[boost::target(vi_out, m_graph)];
      const auto& sp = sv.property<s_prop_t>();
      std::cout << ' ' << sv.get_label();
      if (!inactive) {
        if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
          // in case of a vertex type other than the species or the reaction
          if (sv.get_type() != wcs::Vertex::_species_) continue;
        }
        const auto stoichio = m_graph[vi_out].get_stoichiometry_ratio();
        inactive |= !sp.inc_check(stoichio);
      }
    } // end of for loop over out-edges

    if constexpr (is_bidirectional) {
      std::cout << " from a set of reactants";
      for(const auto vi_in : boost::make_iterator_range(boost::in_edges(vd, m_graph))) {
        const auto& sv = m_graph[boost::source(vi_in, m_graph)];
        const auto& sp = sv.property<s_prop_t>();
        std::cout << ' ' << sv.get_label() << " [" << sp.get_count() << "]";

        if (!inactive) {
          if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
            // in case of a vertex type other than the species or the reaction
            if (sv.get_type() != wcs::Vertex::_species_) continue;
          }
          const auto stoichio = m_graph[vi_in].get_stoichiometry_ratio();
          inactive |= !sp.dec_check(stoichio);
        }
      } // end of for loop over in-edges
    } // end of is_bidirectional

    num_inactive += static_cast<size_t>(inactive);
    if (inactive) {
      m_graph[vd].property<r_prop_t>().set_rate(0.0);
    }

    std::cout << std::endl << "    by the rate " << rp.get_rate()
              << " <= {" << rp.get_rate_formula() << "}" << std::endl;
  } // end of for loop over the reaction list

  std::cout << "Num inactive reactions: "
            << num_inactive << "/" << reaction_list().size() << std::endl;
}

/**@}*/
} // end of namespace wcs
