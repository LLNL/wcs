/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __WCS_BGL_HPP__
#define __WCS_BGL_HPP__

#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
// To suppress the gcc compiler warning 'maybe-uninitialized'
// from the boost graph source code.
// clang does not recognize this particular diagnostic flag.
#include <boost/graph/adjacency_list.hpp>
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#endif

namespace wcs {
/** \addtogroup wcs_global
 *  @{ */

template<typename S = void>
struct adjlist_selector_t {
  using type = ::boost::vecS;
};

template<> struct adjlist_selector_t<::boost::vecS> {
  using type = ::boost::vecS;
};

template<> struct adjlist_selector_t<::boost::listS> {
  using type = ::boost::listS;
};

template<> struct adjlist_selector_t<::boost::setS> {
  using type = ::boost::setS;
};

template<> struct adjlist_selector_t<::boost::multisetS> {
  using type = ::boost::multisetS;
};

template<> struct adjlist_selector_t<::boost::hash_setS> {
  using type = ::boost::hash_setS;
};

#ifndef WCS_VERTEX_LIST_TYPE
using wcs_vertex_list_t = adjlist_selector_t<>::type;
#else
using wcs_vertex_list_t = adjlist_selector_t<WCS_VERTEX_LIST_TYPE>::type;
#endif

#ifndef WCS_OUT_EDGE_LIST_TYPE
using wcs_out_edge_list_t = adjlist_selector_t<>::type;
#else
using wcs_out_edge_list_t = adjlist_selector_t<WCS_OUT_EDGE_LIST_TYPE>::type;
#endif

/**@}*/
} // end of namespace wcs

#endif // __WCS_BGL_HPP__
