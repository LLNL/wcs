/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __WCS_UTILS_GRAPH_FACTORY_HPP__
#define __WCS_UTILS_GRAPH_FACTORY_HPP__

#include <string>
#include <unordered_map>
#include <type_traits>
#include <dlfcn.h> //dlopen
#include "utils/file.hpp"

#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#if defined(WCS_HAS_SBML)
#include <sbml/SBMLTypes.h>
#include <sbml/common/extern.h>
#endif // defined(WCS_HAS_SBML)

// To suppress the gcc compiler warning 'maybe-uninitialized'
// from the boost graph source code.
// clang does not recognize this particular diagnostic flag.
#include <boost/graph/graphml.hpp>
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
#include <bgl.hpp>

#include <reaction_network/vertex_flat.hpp>
#include <reaction_network/vertex.hpp>
#include <reaction_network/edge.hpp>
#include <reaction_network/species.hpp>
#include <reaction_network/reaction.hpp>
#include <utils/detect_methods.hpp>
#include <utils/sbml_utils.hpp>


namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

typedef reaction_rate_t ( * rate_function_pointer)(const std::vector<reaction_rate_t>&);

class GraphFactory {
 public:
  using v_prop_t = wcs::VertexFlat;
  using e_prop_t = wcs::Edge;

  using graph_t = boost::adjacency_list<boost::listS,
                                        boost::listS,
                                        boost::directedS,
                                        v_prop_t,
                                        e_prop_t>;

  using v_desc_t = boost::graph_traits<graph_t>::vertex_descriptor;
  using e_desc_t = boost::graph_traits<graph_t>::edge_descriptor;

  // TODO: gather these JIT-related types into a single file.
  using params_map_t = std::unordered_map <std::string, std::vector<std::string>>;
  using rate_rules_dep_t = std::unordered_map <std::string, std::set<std::string>>;

 public:
  GraphFactory();
  GraphFactory(const GraphFactory &o);
  const graph_t &graph() const;
  bool read_graphml(const std::string &ifn);
  /** Export the internal adjacency list to g, which might be of a different
      type in terms of the random accessibility but of a compatible one (G). */
  template<typename G> void copy_to(G& g) const;

  #if defined(WCS_HAS_SBML)
  template<typename G>
  void convert_to(const LIBSBML_CPP_NAMESPACE::Model& model, G& g,
    const std::string& library_file,
    const params_map_t& dep_params_f,
    const params_map_t& dep_params_nf,
    const rate_rules_dep_t& rate_rules_dep_map) const;
  #endif // defined(WCS_HAS_SBML)

  template<typename G>
  typename std::shared_ptr<G> make_graph() const;

 private:
  void setup_dynamic_property ();
  graph_t m_g;
  boost::dynamic_properties m_dp;
};


template<typename G> void GraphFactory::copy_to(G& g) const
{
  using v_new_desc_t = typename boost::graph_traits<G>::vertex_descriptor;
  std::unordered_map<v_desc_t, v_new_desc_t> v_desc_map;

  if constexpr (has_reserve_for_vertex_list<G>::value) {
    g.m_vertices.reserve(m_g.m_vertices.size());
  }

  using directed_category = typename boost::graph_traits<G>::directed_category;
  constexpr bool is_bidirectional
    = std::is_same<directed_category, boost::bidirectional_tag>::value;

  typename boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;

  using in_reactions = std::unordered_map<v_desc_t, unsigned int> ;
  typename in_reactions::const_iterator idit;
  in_reactions in_degree_reaction;

  //calculate the in_edges of each reaction
  for (boost::tie(vi, vi_end) = boost::vertices(m_g); vi != vi_end; ++vi) {
    const v_prop_t& v = m_g[*vi];
    const auto vt = static_cast<v_prop_t::vertex_type>(v.get_typeid());

    if (vt == v_prop_t::_species_) {
      for (const auto ei_out :
           boost::make_iterator_range(boost::out_edges(*vi, m_g)))
      {
        const auto vd_r = boost::target(ei_out, m_g);
        in_degree_reaction[vd_r] += 1;
      }
    }
  }


  for (boost::tie(vi, vi_end) = boost::vertices(m_g); vi != vi_end; ++vi) {
    const v_prop_t& v = m_g[*vi];
    const auto vt = static_cast<v_prop_t::vertex_type>(v.get_typeid());
    const auto num_out_edges = boost::out_degree(*vi, m_g);
    v_new_desc_t vd;

    if (vt == v_prop_t::_species_) {
      vd = boost::add_vertex(wcs::Vertex{v, g}, g);
      if constexpr (has_reserve_for_vertex_out_edges<G>::value) {
        g.m_vertices[vd].m_out_edges.reserve(num_out_edges);
      } else { (void) num_out_edges; }

      if constexpr (has_reserve_for_vertex_in_edges<G>::value &&
                    is_bidirectional)
      {
        g.m_vertices[vd].m_in_edges.reserve(num_in_edges_to_reserve);
      }
      v_desc_map[*vi] = vd;
    } else if (vt == v_prop_t::_reaction_) {
      idit = in_degree_reaction.find(*vi) ;
      if (idit != in_degree_reaction.end() || num_out_edges != 0u) {
        vd = boost::add_vertex(wcs::Vertex{v, g}, g);
        if constexpr (has_reserve_for_vertex_out_edges<G>::value) {
          g.m_vertices[vd].m_out_edges.reserve(num_out_edges);
        }
        if constexpr (has_reserve_for_vertex_in_edges<G>::value &&
                      is_bidirectional)
        {
          g.m_vertices[vd].m_in_edges.reserve(num_in_edges_to_reserve);
        }
        v_desc_map[*vi] = vd;
      }
    }
  }

  typename boost::graph_traits<graph_t>::edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = boost::edges(m_g); ei != ei_end; ++ei) {
    const v_desc_t u = boost::source(*ei, m_g),
                   v = boost::target(*ei, m_g);

    v_new_desc_t u_new = v_desc_map.at(u);
    v_new_desc_t v_new = v_desc_map.at(v);
    const e_prop_t& e = m_g[*ei];
    boost::add_edge(u_new, v_new, e_prop_t{e}, g);
  }
}

#if defined(WCS_HAS_SBML)
/// Create a Boost graph out of an SBML model
template<typename G> void
GraphFactory::convert_to(
  const LIBSBML_CPP_NAMESPACE::Model& model, G& g,
  const std::string& library_file,
  const params_map_t& dep_params_f,
  const params_map_t& dep_params_nf,
  const rate_rules_dep_t& rate_rules_dep_map) const
{
  using v_new_desc_t = typename boost::graph_traits<G>::vertex_descriptor;
  using e_new_desc_t = typename boost::graph_traits<G>::edge_descriptor;


  if constexpr (has_reserve_for_vertex_list<G>::value) {
    g.m_vertices.reserve(m_g.m_vertices.size());
  }

  using directed_category = typename boost::graph_traits<G>::directed_category;
  constexpr bool is_bidirectional
    = std::is_same<directed_category, boost::bidirectional_tag>::value;

  const LIBSBML_CPP_NAMESPACE::ListOfReactions* reaction_list
    = model.getListOfReactions();

  if  (reaction_list == nullptr) {
    WCS_THROW("Invalid reaction_list pointer");
    return;
  }
  unsigned int num_reactions = reaction_list->size();
  const LIBSBML_CPP_NAMESPACE::ListOfSpecies* species_list
    = model.getListOfSpecies();

  using species_added_t = std::unordered_map<std::string, v_new_desc_t>;
  typename species_added_t::const_iterator it;
  species_added_t smap;

  using edges_added_t = std::unordered_map<std::string, e_new_desc_t>;
  typename edges_added_t::const_iterator eit;
  edges_added_t emap;

  using undeclared_reactants_t = std::unordered_set<std::string>;
  typename undeclared_reactants_t::const_iterator urit;
  undeclared_reactants_t urset;

  using  all_species_t = std::unordered_set<std::string>;
  typename all_species_t::const_iterator aspit;
  all_species_t aspset;

  // Create an unordered_set for all model species
  unsigned int num_species = species_list->size();
  for (unsigned int si = 0u; si < num_species; si++) {
    aspset.insert(species_list->get(si)->getIdAttribute());
  }

  std::string library_name = library_file;
  #if !defined(WCS_HAS_EXPRTK)
  if (library_name.find_first_of("/") == std::string::npos) {
    // Only if filename contains "/", then it is interpreted as a (relative or
    // absolute) pathname. Otherwise, dlopen searches the library in various
    // paths listed in  DT_PATH, LD_LIBRARY_PATH, DT_RUNPATH and system
    // directories
    library_name = "./" + library_name;
  }
  // Implement the interface to the dynamic linking loader
  void* handle = dlopen(library_name.c_str(), RTLD_LAZY);

  if (!handle) {
    WCS_THROW("Cannot open library '" + library_name + "': " \
              + std::string(dlerror()) + "\n");
    return;
  }

  // reset errors
  dlerror();
  #endif // !defined(WCS_HAS_EXPRTK)

  // Add reactions
  for (unsigned int ri = 0u; ri < num_reactions; ri++) {
    std::unordered_set<std::string> reaction_in_species;
    std::unordered_set<std::string>::const_iterator risit;
    if (reaction_list->get(ri) == nullptr) {
      WCS_THROW("Invalid reaction pointer");
      return;
    }
    const auto &reaction = *(reaction_list->get(ri));

    #if !defined(WCS_HAS_EXPRTK)
    std::string rate_function = std::string("wcs__rate_" + reaction.getIdAttribute());

    void * const reaction_function_rate = dlsym(handle, rate_function.c_str());
    const char *dlsym_error = dlerror();
    if (dlsym_error) {
      WCS_THROW("Cannot load symbol " + reaction.getIdAttribute() + ": " + dlsym_error + "\n");
      dlclose(handle);
      return;
    }

    std::function<reaction_rate_t(const std::vector<reaction_rate_t>&)> f_rate =
      reinterpret_cast<rate_function_pointer>(reaction_function_rate);
    wcs::Vertex v(model, reaction, g, f_rate);
    #else
    wcs::Vertex v(model, reaction, g, NULL);
    #endif // !defined(WCS_HAS_EXPRTK)

    v_new_desc_t vd = boost::add_vertex(v, g);
    unsigned int num_reactants = reaction.getNumReactants();
    unsigned int num_products = reaction.getNumProducts();
    unsigned int num_modifiers = reaction.getNumModifiers();

    if constexpr (has_reserve_for_vertex_out_edges<G>::value) {
        g.m_vertices[vd].m_out_edges.reserve(num_products);
    }
    if constexpr (has_reserve_for_vertex_in_edges<G>::value &&
                  is_bidirectional)
    {
      g.m_vertices[vd].m_in_edges.reserve(num_reactants);
    }
    // Add reactants species
    for (unsigned int si = 0u; si < num_reactants; si++) {
      const auto &reactant = *(reaction.getReactant(si));

      std::string s_label
        = species_list->get(reactant.getSpecies())->getIdAttribute();
      it = smap.find(s_label) ;
      v_new_desc_t vds;
      reaction_in_species.insert(s_label);

      if (it == smap.cend()) {
        wcs::Vertex vs(model, *species_list->get(reactant.getSpecies()), g);
        //WCS_THROW("TESTT");
        vds = boost::add_vertex(vs, g);
        smap.insert(std::make_pair(s_label,vds));
      } else {
        vds = it->second;
      }

      std::string e_label = g[vds].get_label() + '|' + g[vd].get_label();
      if (!species_list->get(reactant.getSpecies())->getConstant()) { //Add as a reactant
        eit = emap.find(e_label);

        if (eit == emap.cend()) {
          const auto ret = boost::add_edge(vds, vd, g);

          if (!ret.second) {
            WCS_THROW("Please check the reactions in your SBML file");
            return;
          }
          g[ret.first].set_stoichiometry_ratio(reactant.getStoichiometry());
          g[ret.first].set_label(e_label);
          emap.insert(std::make_pair(e_label, ret.first));
        } else {
          const auto& edge_found = eit->second;
          stoic_t new_stoich = g[edge_found].get_stoichiometry_ratio()
                            + reactant.getStoichiometry();
          g[edge_found].set_stoichiometry_ratio(new_stoich);
        }
      } else {  //Add as a modifier
        const auto ret = boost::add_edge(vds, vd, g);

        if (!ret.second) {
          WCS_THROW("Please check the reactions in your SBML file");
          return;
        }
        g[ret.first].set_stoichiometry_ratio(0);
        g[ret.first].set_label(e_label);
        emap.insert(std::make_pair(e_label, ret.first));
      }  
    }

    // Add modifiers species
    for (unsigned int si = 0u; si < num_modifiers; si++) {
      const auto &modifier = *(reaction.getModifier(si));

      std::string s_label
        = species_list->get(modifier.getSpecies())->getIdAttribute();
      it = smap.find(s_label) ;
      v_new_desc_t vds;
      reaction_in_species.insert(s_label);

      if (it == smap.cend()) {
        wcs::Vertex vs(model, *species_list->get(modifier.getSpecies()), g);
        vds = boost::add_vertex(vs, g);
        smap.insert(std::make_pair(s_label,vds));
      } else {
        vds = it->second;
      }

      std::string e_label = g[vds].get_label() + '|' + g[vd].get_label();
      const auto ret = boost::add_edge(vds, vd, g);

      if (!ret.second) {
        WCS_THROW("Please check the reactions in your SBML file");
        return;
      }
      g[ret.first].set_stoichiometry_ratio(0);
      g[ret.first].set_label(e_label);
      emap.insert(std::make_pair(e_label, ret.first));
    }

    sbml_utils sbml_o;
    urset = sbml_o.find_undeclared_species_in_reaction_formula(model, reaction);

    // Check if all the elements of the undeclared elements are actually species
    for  (const std::string& x: urset) {
      std::string s_label = x;
      aspit = aspset.find(s_label);
      if (aspit == aspset.cend()) {
        WCS_THROW("Unknown element " + s_label + " in the reaction " +
                  reaction.getIdAttribute() + " of your SBML file");
        return;
      }
    }

    // Add undeclared reactants species in the rate formula
    for  (const std::string& x: urset) {
      std::string s_label = x;
      it = smap.find(s_label) ;
      v_new_desc_t vds;
      reaction_in_species.insert(s_label);

      if (it == smap.cend()) {
        wcs::Vertex vs(model, *species_list->get(x), g);
        vds = boost::add_vertex(vs, g);
        smap.insert(std::make_pair(s_label,vds));
      } else {
        vds = it->second;
      }

      std::string e_label = g[vds].get_label() + '|' + g[vd].get_label();
      const auto ret = boost::add_edge(vds, vd, g);

      if (!ret.second) {
        WCS_THROW("Please check the reactions in your SBML file");
        return;
      }
      g[ret.first].set_stoichiometry_ratio(0);
      g[ret.first].set_label(e_label);
      emap.insert(std::make_pair(e_label, ret.first));
    }

    // Add transient parameters
    typename params_map_t::const_iterator pit_f, pit_nf;
    typename rate_rules_dep_t::const_iterator rrdit;
    //transient parameters in formula
    pit_f = dep_params_f.find(reaction.getIdAttribute());
    if (pit_f != dep_params_f.cend()) {
      const std::vector<std::string>& params_f = pit_f->second;
      std::vector<std::string>::const_iterator itf;
      for (itf = params_f.cbegin(); itf != params_f.cend(); itf++){
        rrdit = rate_rules_dep_map.find(*itf);
        if (rrdit != rate_rules_dep_map.cend()) {
          const std::set<std::string>& params_rr = rrdit->second;
          for (const auto& x: params_rr) {
            it = smap.find(x);
            v_new_desc_t vds;
            risit = reaction_in_species.find(x);
            if (risit == reaction_in_species.cend()) {
              reaction_in_species.insert(x);
              if (it == smap.cend()) {
                wcs::Vertex vs(model, *species_list->get(x), g);
                vds = boost::add_vertex(vs, g);
                smap.insert(std::make_pair(x,vds));
              } else {
                vds = it->second;
              }

              std::string e_label = g[vds].get_label() + '|' + g[vd].get_label();
              const auto ret = boost::add_edge(vds, vd, g);

              if (!ret.second) {
                WCS_THROW("Please check the reactions in your SBML file");
                return;
              }
              g[ret.first].set_stoichiometry_ratio(0);
              g[ret.first].set_label(e_label);
              emap.insert(std::make_pair(e_label, ret.first));
            }
          }
        }
      }
    }

    //transient parameters not in formula
    pit_nf = dep_params_nf.find(reaction.getIdAttribute());
    if (pit_nf != dep_params_nf.cend()) {
      const std::vector<std::string>& params_nf = pit_nf->second;
      std::vector<std::string>::const_iterator itf;
      for (itf = params_nf.cbegin(); itf != params_nf.cend(); itf++){
        rrdit = rate_rules_dep_map.find(*itf);
        if (rrdit != rate_rules_dep_map.cend()) {
          const std::set<std::string>& params_rr = rrdit->second;
          for (const auto& x: params_rr) {
            it = smap.find(x);
            v_new_desc_t vds;
            risit = reaction_in_species.find(x);
            if (risit == reaction_in_species.cend()) {
              if (it == smap.cend()) {
                wcs::Vertex vs(model, *species_list->get(x), g);
                vds = boost::add_vertex(vs, g);
                smap.insert(std::make_pair(x,vds));
              } else {
                vds = it->second;
              }

              reaction_in_species.insert(x);
              std::string e_label = g[vds].get_label() + '|' + g[vd].get_label();
              const auto ret = boost::add_edge(vds, vd, g);

              if (!ret.second) {
                WCS_THROW("Please check the reactions in your SBML file");
                return;
              }
              g[ret.first].set_stoichiometry_ratio(0);
              g[ret.first].set_label(e_label);
              emap.insert(std::make_pair(e_label, ret.first));
            }
          }
        }
      }
    }

    // Add products species
    for (unsigned int si = 0u; si < num_products; si++) {
      const auto &product = *(reaction.getProduct(si));

      std::string s_label
        = species_list->get(product.getSpecies())->getIdAttribute();
      it = smap.find(s_label);
      v_new_desc_t vds;

      if (it == smap.cend()) {
        wcs::Vertex vs(model, *species_list->get(product.getSpecies()), g);
        vds = boost::add_vertex(vs, g);
        smap.insert(std::make_pair(s_label,vds));
      } else {
        vds = it->second;
      }

      std::string e_label = g[vd].get_label() + '|' + g[vds].get_label();

      if (!species_list->get(product.getSpecies())->getConstant()) { //Add as a product
        eit = emap.find(e_label);
        if (eit == emap.cend()) {
          const auto ret = boost::add_edge(vd, vds, g);

          if (!ret.second) {
            WCS_THROW("Please check the reactions in your SBML file");
            return;
          }
          g[ret.first].set_stoichiometry_ratio(product.getStoichiometry());
          g[ret.first].set_label( g[vd].get_label() + '|' + g[vds].get_label());
          emap.insert(std::make_pair(e_label, ret.first));
        } else {
          const auto& edge_found = eit->second;
          stoic_t new_stoich = g[edge_found].get_stoichiometry_ratio()
                            + product.getStoichiometry();
          g[edge_found].set_stoichiometry_ratio(new_stoich);
        }
      } else { //Add as a modifier
        const auto ret = boost::add_edge(vds, vd, g);

        if (!ret.second) {
          WCS_THROW("Please check the reactions in your SBML file");
          return;
        }
        g[ret.first].set_stoichiometry_ratio(0);
        g[ret.first].set_label(e_label);
        emap.insert(std::make_pair(e_label, ret.first));
      }
    }
  }
}
#endif // defined(WCS_HAS_SBML)


/**
 * Load the graph from GraphML file and return the shared pointer of the
 * graph object of type G.
 */
template<typename G>
typename std::shared_ptr<G> GraphFactory::make_graph() const
{
  std::shared_ptr<G> g = std::make_shared<G>();
  copy_to(*g);
  return g;
}

/**@}*/
} // end of namespace wcs
#endif // __WCS_UTILS_GRAPH_FACTORY_HPP__
