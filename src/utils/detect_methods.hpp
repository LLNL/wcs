#ifndef __WCS_UTILS_DETECT_METHODS_HPP__
#define __WCS_UTILS_DETECT_METHODS_HPP__

#include <type_traits>
#include <boost/graph/adjacency_list.hpp>

namespace wcs {
/** \addtogroup wcs_utils
 *  *  @{ */

/**
 * Check if the vertex list of the graph of type G indicated by the first
 * template parameter has the interface `reserve()` to pre-allocate the
 * space. The check is done at compile time incurring no run-time overhead.
 */
template <typename Graph, typename Size = size_t>
struct has_reserve_for_vertex_list {
  template <typename G, typename S>
  static constexpr
  decltype(std::declval<G>().m_vertices.reserve(std::declval<S>()), bool())
  test_reserve(int) {
    return true;
  }

  template <typename G, typename S>
  static constexpr bool test_reserve(...) {
    return false;
  }

  static constexpr bool value = test_reserve<Graph, Size>(int());
};

/**
 * Check if the edge list of the graph of type G indicated by the first
 * template parameter has the interface `reserve()` to pre-allocate the
 * space. The check is done at compile time incurring no run-time overhead.
 */
template <typename Graph, typename Size = size_t>
struct has_reserve_for_edge_list {
  template <typename G, typename S>
  static constexpr
  decltype(std::declval<G>().m_edges.reserve(std::declval<S>()), bool())
  test_reserve(int) {
    return true;
  }

  template <typename G, typename S>
  static constexpr bool test_reserve(...) {
    return false;
  }

  static constexpr bool value = test_reserve<Graph, Size>(int());
};

/**
 * Check if the out-edge list of the vertex in the vertex list of the graph
 * of type G indicated by the first template parameter has the interface
 * `reserve()` to pre-allocate the space. The check is done at compile time
 * incurring no run-time overhead.
 */
template <typename Graph, typename Size = size_t>
struct has_reserve_for_vertex_out_edges {
  template <typename G, typename S>
  static constexpr
  decltype(std::declval<G>().m_vertices.back().m_out_edges.reserve(std::declval<S>()), bool())
  test_reserve(int) {
    return true;
  }

  template <typename G, typename S>
  static constexpr bool test_reserve(...) {
    return false;
  }

  static constexpr bool value = test_reserve<Graph, Size>(int());
};

/**
 * Check if the in-edge list of the vertex in the vertex list of the graph
 * of type G indicated by the first template parameter has the interface
 * `reserve()` to pre-allocate the space. The check is done at compile time
 * incurring no run-time overhead.
 */
template <typename Graph, typename Size = size_t>
struct has_reserve_for_vertex_in_edges {
  template <typename G, typename S>
  static constexpr
  decltype(std::declval<G>().m_vertices.back().m_in_edges.reserve(std::declval<S>()), bool())
  test_reserve(int) {
    return true;
  }

  template <typename G, typename S>
  static constexpr bool test_reserve(...) {
    return false;
  }

  static constexpr bool value = test_reserve<Graph, Size>(int());
};
/**@}*/
} // end of namespace wcs
#endif // __WCS_UTILS_DETECT_METHODS_HPP__
