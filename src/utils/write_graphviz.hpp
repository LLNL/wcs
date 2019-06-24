#ifndef __WCS_UTILS_WRITE_GRAPHVIZ_HPP__
#define __WCS_UTILS_WRITE_GRAPHVIZ_HPP__
#include <boost/graph/graphviz.hpp>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <algorithm>
#include <type_traits>
#include <unordered_map>


namespace wcs {
/** \addtogroup wcs_utils
 *  *  @{ */

/**
 * Check if the associative container has the interface `at()` to access the value
 * by a key. The check is done at compile time incurring no run-time overhead.
 */
template <typename Map, typename Key>
struct has_at {
  template <typename M, typename K>
  static constexpr
  decltype(std::declval<M>().at(std::declval<K>()), bool())
  test_at(int) {
    return true;
  }

  template <typename M, typename K>
  static constexpr bool test_at(...) {
    return false;
  }

  static constexpr bool value = test_at<Map, Key>(int());
};


/**
 * Write the graph out to the output stream in the graphviz format.
 */
template <typename G, typename VIdxMap>
std::ostream& write_graphviz(std::ostream& os, const G& g, const VIdxMap& v_idx_map)
{
  using std::operator<<;

  using directed_category = typename boost::graph_traits<G>::directed_category;
  constexpr bool directed = std::is_same<directed_category, boost::directed_tag>::value
                         || std::is_same<directed_category, boost::bidirectional_tag>::value;

  os << (directed? "digraph G {" : "graph G {") << std::endl
     << "node [margin=0 fontcolor=blue fontsize=32 width=0.5 shape=circle style=filled]"
     << std::endl;

  using v_desc_t = typename boost::graph_traits<G>::vertex_descriptor;

  std::function<size_t (const VIdxMap&, const v_desc_t&)> get_v_index
    = [] (const VIdxMap& vidx_map, const v_desc_t& v) {
      if constexpr (has_at<VIdxMap, v_desc_t>::value) {
        return vidx_map.at(v);
      } else {
        return vidx_map[v];
      }
    };

  { // print vertices
    using v_prop_t = typename boost::vertex_bundle_type<G>::type;
    //using v_prop_t = typename std::remove_cv<typename std::remove_reference<decltype(g[*vi])>::type>::type;
    //using v_prop_t = typename G::vertex_property_type;

    // Use detection idiom to see if at() is available.
    // If so, replace g[v] with g.at(v), and this allows const 
    // v_idx_map with std containers, which does not define
    // const access by the [] operator.
    std::function<void (const v_desc_t&)> vertex_writer
      = [&] (const v_desc_t& v) {
        static std::map<typename v_prop_t::vertex_type, std::string> colormap {
          {v_prop_t::_undefined_, "white"},
          {v_prop_t::_species_,   "lightsteelblue"},
          {v_prop_t::_reaction_,  "plum1"}
        };
  
        os << "  " << get_v_index(v_idx_map, v)
           << " [label=\"" << boost::escape_dot_string(g[v].get_label())
           << "\" fillcolor=\"" << colormap.at(g[v].get_type())
           << "\"];" << std::endl;
      };
  
    typename boost::graph_traits<G>::vertex_iterator vi, vi_end;
    boost::tie(vi, vi_end) = boost::vertices(g);

    std::for_each(vi, vi_end, vertex_writer);
  }

  { // print edges
    using e_desc_t = typename boost::graph_traits<G>::edge_descriptor;
  
    constexpr const char* edge_delimiter = directed? "->" : "--";

    std::function<void (const e_desc_t&)> edge_writer
      = [&] (const e_desc_t& e) {
        os << "  " << get_v_index(v_idx_map, boost::source(e, g))
           << edge_delimiter
           << get_v_index(v_idx_map, boost::target(e, g))
           << " [taillabel=\"" << g[e].get_stoichiometry_ratio() << "\"];"
           << std::endl;
      };
  
    typename boost::graph_traits<G>::edge_iterator ei, ei_end;
    boost::tie(ei, ei_end) = boost::edges(g);

    std::for_each(ei, ei_end, edge_writer);
  }

  os << "}" << std::endl;

  return os;
}


template <typename G>
std::ostream& write_graphviz(std::ostream& os, const G& g)
{
  using rand_access
    = typename boost::detail::is_random_access<typename G::vertex_list_selector>::type;
   
  if constexpr (rand_access::value) {
    using v_index_map_t
      = typename boost::property_map<G, boost::vertex_index_t>::const_type;
    v_index_map_t vidx_map = boost::get(boost::vertex_index, g);
    ::wcs::write_graphviz(os, g, vidx_map);
  } else {
    using v_desc_t = typename boost::graph_traits<G>::vertex_descriptor;
    std::unordered_map<v_desc_t, size_t> v_idx_map;
    for (auto u : boost::make_iterator_range(boost::vertices(g))) {
      v_idx_map[u] = v_idx_map.size();
    }
    ::wcs::write_graphviz(os, g, v_idx_map);
  }
  return os;
}


template <typename G>
std::ostream& operator<<(std::ostream& os, const G& g)
{
  ::wcs::write_graphviz(os, g);
  return os;
}


template <typename G>
bool write_graphviz(const std::string& out_filename, const G& g)
{
  bool ok = true;
  std::ofstream ofs (out_filename.c_str());

  try {
    ::wcs::write_graphviz(ofs, g);
  } catch (boost::graph_exception &e) {
    using std::operator<<;
    std::cerr << e.what () << std::endl;
    ok = false;
  }
  ofs.close ();
  return ok;
}

/**@}*/
} // end of namespace wcs

#endif // __WCS_UTILS_WRITE_GRAPHVIZ_HPP__
