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

#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
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

#include <sbml/SBMLTypes.h>   ///Konstantia added
#include <sbml/common/extern.h>


namespace wcs {
/** \addtogroup wcs_utils
 *  *  @{ */

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

 public:
  GraphFactory();
  GraphFactory(const GraphFactory &o);
  const graph_t &graph() const;
  bool read_graphml(const std::string &ifn);
  /** Export the internal adjacency list to g, which might be of a different
      type in terms of the random accessibility but of a compatible one (G). */
  template<typename G> void copy_to(G& g) const;
  template<typename G> void convert_to(libsbml::Model& model, G& g) const; ///Konstantia
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
      if constexpr (has_reserve_for_vertex_in_edges<G>::value && is_bidirectional) {
        g.m_vertices[vd].m_in_edges.reserve(num_in_edges_to_reserve);
      }
    } else {
      vd = boost::add_vertex(wcs::Vertex{v, g}, g);
      if constexpr (has_reserve_for_vertex_out_edges<G>::value) {
        g.m_vertices[vd].m_out_edges.reserve(num_out_edges);
      }
      if constexpr (has_reserve_for_vertex_in_edges<G>::value && is_bidirectional) {
        g.m_vertices[vd].m_in_edges.reserve(num_in_edges_to_reserve);
      }
    }
    v_desc_map[*vi] = vd;
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

/// Convert SBML output into a Boost graph 
template<typename G> void GraphFactory::convert_to(libsbml::Model& model, G& g) const
{
  using v_new_desc_t = typename boost::graph_traits<G>::vertex_descriptor;

  if constexpr (has_reserve_for_vertex_list<G>::value) {
    g.m_vertices.reserve(m_g.m_vertices.size());
  }

  using directed_category = typename boost::graph_traits<G>::directed_category;
  constexpr bool is_bidirectional
    = std::is_same<directed_category, boost::bidirectional_tag>::value;


  libsbml::ListOfReactions* reactionslist = model.getListOfReactions();

  if  (reactionslist == nullptr) {
    WCS_THROW("Invalid reactionslist pointer");
    return; 
  }
  unsigned int reactionsSize = reactionslist->size();
  libsbml::ListOfSpecies* specieslist = model.getListOfSpecies();

  using species_added = std::unordered_map<std::string, v_new_desc_t>;
  typename species_added::const_iterator it;
  species_added smap;

  using edges_added = std::unordered_map<std::string, double>;
  edges_added::const_iterator eit;
  edges_added emap;

  /// Add reactions
  for (unsigned int ri=0u; ri < reactionsSize; ri++) {
    if (reactionslist->get(ri) == nullptr) {
      WCS_THROW("Invalid reaction pointer");
      return;
    } 
    const auto &reaction = *(reactionslist->get(ri));
        
    wcs::Vertex v(model, reaction, g);
    
    v_new_desc_t vd = boost::add_vertex(v, g);
    unsigned int reactantsSize = reaction.getNumReactants();
    unsigned int productsSize = reaction.getNumProducts();
    if constexpr (has_reserve_for_vertex_out_edges<G>::value) {
        g.m_vertices[vd].m_out_edges.reserve(productsSize);
    }
    if constexpr (has_reserve_for_vertex_in_edges<G>::value && is_bidirectional) {
      g.m_vertices[vd].m_in_edges.reserve(reactantsSize);
    }

    /// Add reactants species
    for (unsigned int si=0u; si < reactantsSize; si++) {
      const auto &reactant = *(reaction.getReactant(si));
      
      std::string s_label = specieslist->get(reactant.getSpecies())->getIdAttribute();
      it = smap.find(s_label) ;
      v_new_desc_t vds;
      if (it == smap.end()) {
        wcs::Vertex vs(*specieslist->get(reactant.getSpecies()), g); 
        vds = boost::add_vertex(vs, g); 
        smap.insert(std::make_pair(s_label,vds));
      } else {
        vds = it->second;
      }

      std::string e_label = g[vds].get_label() + g[vd].get_label();
      eit = emap.find(e_label);
      if (eit == emap.end()) {
        const auto ret = boost::add_edge(vds, vd, g);
      
        if (!ret.second) { 
          WCS_THROW("Please check the reactions in your SBML file");
          return;
        }
        g[ret.first].set_stoichiometry_ratio(reactant.getStoichiometry()); 
        g[ret.first].set_label(e_label);
        emap.insert(std::make_pair(e_label,reactant.getStoichiometry()));
      }  
    }
    
    /// Add products species
    for (unsigned int si=0u; si < productsSize; si++) {
      const auto &product = *(reaction.getProduct(si));
   
      std::string s_label = specieslist->get(product.getSpecies())->getIdAttribute();
      it = smap.find(s_label) ;
      v_new_desc_t vds;
      if (it == smap.end()) {
        wcs::Vertex vs(*specieslist->get(product.getSpecies()), g); 
        vds = boost::add_vertex(vs, g); 
        smap.insert(std::make_pair(s_label,vds));
      } else {
        vds = it->second;
      }

      std::string e_label = g[vd].get_label() + g[vds].get_label();
      eit = emap.find(e_label);
      if (eit == emap.end()) {
        const auto ret = boost::add_edge(vd, vds, g);
        
        if (!ret.second) {
          WCS_THROW("Please check the reactions in your SBML file");
          return;
        }
        g[ret.first].set_stoichiometry_ratio(product.getStoichiometry()); 
        g[ret.first].set_label( g[vd].get_label() + g[vds].get_label()); 
        emap.insert(std::make_pair(e_label,product.getStoichiometry()));
      }
    }
    
  }  
}


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
