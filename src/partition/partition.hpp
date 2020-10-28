/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __WCS_PARTITION_PARTITION_HPP__
#define __WCS_PARTITION_PARTITION_HPP__

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

#include <reaction_network/network.hpp>
#include <reaction_network/edge_weighted.hpp>

namespace wcs {
/** \addtogroup wcs_partition
 *  @{ */

class Partition
{
 public:
  /** The type of the BGL graph to represent an intermediate graph transformed
   *  from a reaction graph for partitioning.
   */
  using intm_graph_t = boost::adjacency_list<
                       wcs_out_edge_list_t,
                       wcs_vertex_list_t,
                       boost::bidirectionalS,
                       wcs::Vertex, // vertex property bundle
                       wcs::Edge_Weighted, // edge property bundle
                       boost::no_property,
                       boost::vecS>;

  using vd_t = typename boost::graph_traits<intm_graph_t>::vertex_descriptor;
  using pid_t = wcs::partition_id_t;


  /**
   * Convert a reaction graph to an intermediate one to faciliate partitioning.
   * The graph before conversion is an unweighted bipartite graph that consists
   * of reaction vertices as well as non-reaction vertices such as species. The
   * new graph is weighted and not bipartite. It does not have a reaction
   * vertex. Instead, a set of edges from each input of a reaction to each
   * output of it substitutes the reaction vertex.
   */
  void make_intermediate_graph(const wcs::Network& net);

  //std::vector<pid_t> find_vertex_assignment(vd_t vd) const;

  const intm_graph_t& intm_graph() const;
  intm_graph_t& intm_graph();

 protected:
  intm_graph_t m_g_intm;
};


/**@}*/
} // end of namespace wcs
#endif // __WCS_PARTITION_PARTITION_HPP__
