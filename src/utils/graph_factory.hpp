#ifndef __WCS_UTILS_GRAPH_FACTORY_HPP__
#define __WCS_UTILS_GRAPH_FACTORY_HPP__

#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
// To suppress the gcc compiler warning 'maybe-uninitialized'
// from the boost graph source code.
// clang does not recognize this particular diagnostic flag.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphml.hpp>
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#include <reaction_network/vertex_flat.hpp>
#include <reaction_network/vertex.hpp>
#include <reaction_network/edge.hpp>
#include <reaction_network/species.hpp>
#include <reaction_network/reaction.hpp>
#include <utils/detect_methods.hpp>
#include <string>
#include <unordered_map>
#include <type_traits>


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
  template<typename G> void copy_to(G& g);

 private:
  void setup_dynamic_property ();
  graph_t m_g;
  boost::dynamic_properties m_dp;
};


template<typename G> void GraphFactory::copy_to(G& g)
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
    //std::cerr << "adding edge from "
    //          << g[u_new].get_label() << " to "
    //          << g[v_new].get_label() << std::endl;
  }
}

/**@}*/
} // end of namespace wcs
#endif // __WCS_UTILS_GRAPH_FACTORY_HPP__
