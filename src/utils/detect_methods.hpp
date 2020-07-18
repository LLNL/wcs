/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

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
struct has_reserve_for_vertex_list
{
 private:
  typedef std::true_type yes;
  typedef std::false_type no;

  template <typename G, typename S>
  static auto test(int) ->
    decltype(std::declval<G>().m_vertices.reserve(std::declval<S>()), yes());

  template <typename G, typename S>
  static no test(...);

 public:
  using type = decltype(test<Graph, Size>(0));
  static constexpr bool value = std::is_same<type, yes>::value;
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

/**
 * Check if the associative container has the interface `at()` to access the value
 * by a key. The check is done at compile time incurring no run-time overhead.
 */
template <typename Map, typename Key>
struct has_at {
private:
  typedef std::true_type yes;
  typedef std::false_type no;

  template <typename M, typename K>
  static auto test(int) ->
    decltype(std::declval<M>().at(std::declval<K>()), yes());

  template <typename M, typename K>
  static no test(...);

public:
  using type = decltype(test<Map,Key>(0));
  static constexpr bool value = std::is_same<type, yes>::value;
};

/// Detect if edge type has a member get_weight()
template <typename E>
struct has_get_weight {
private:
  typedef std::true_type yes;
  typedef std::false_type no;

  template <typename EE>
  static auto test(int) ->
    decltype(std::declval<EE>().get_weight() ==
             std::declval<reaction_rate_t>(), yes());

  template <typename EE>
  static no test(...);

public:
  using type = decltype(test<E>(0));
  static constexpr bool value = std::is_same<type, yes>::value;
};

/// Detect if edge type has a member set_weight()
template <typename E>
struct has_set_weight {
private:
  typedef std::true_type yes;
  typedef std::false_type no;

  template <typename EE>
  static auto test(int) ->
    decltype(std::declval<EE>().set_weight(std::declval<reaction_rate_t>()),
             yes());

  template <typename EE>
  static no test(...);

public:
  using type = decltype(test<E>(0));
  static constexpr bool value = std::is_same<type, yes>::value;
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_UTILS_DETECT_METHODS_HPP__
