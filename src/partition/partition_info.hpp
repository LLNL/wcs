/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __WCS_PARTITION_PARTITION_INFO_HPP__
#define __WCS_PARTITION_PARTITION_INFO_HPP__

#include <bgl.hpp>

#include <reaction_network/network.hpp>
#include <reaction_network/edge_weighted.hpp>

namespace wcs {
/** \addtogroup wcs_partition
 *  @{ */

class Partition_Info
{
 public:
  using graph_t = wcs::Network::graph_t;
  using v_desc_t = wcs::Network::v_desc_t;
  using v_prop_t = wcs::Network::v_prop_t;
  using v_iter_t = wcs::Network::v_iter_t;
  using e_desc_t = wcs::Network::e_desc_t;
  using e_prop_t = wcs::Network::e_prop_t;
  using e_iter_t = wcs::Network::e_iter_t;
  using r_prop_t = wcs::Network::r_prop_t;

  using directed_category = boost::graph_traits<graph_t>::directed_category;
  static constexpr bool is_bidirectional
    = std::is_same<directed_category, boost::bidirectional_tag>::value;
  static constexpr auto n_types = static_cast<size_t>(v_prop_t::_num_vertex_types_ - 1);

  /// Vertex, degree
  using degree_t = std::pair<v_desc_t, size_t>; // TODO: need to be a set
  /// Type for the container of the maximum degree per partition
  using p_max_t = std::vector<degree_t>;
  /// Type for the container of the total degree per partition
  using p_sum_t = std::vector<size_t>;
  /// Type for the computational load
  using load_t = reaction_rate_t;
  /// Type for the container of the compute loads per partition
  using p_load_t = std::vector<load_t>;

  /// Sum of the weight of edges connecting vertices across partitions
  using comm_volume_t = Edge_Weighted::edge_weight_t;
  //using comm_volume_t = stoic_t;
  /// Out-going communication volume to neighbor partitions
  using part_edges_t = std::map<partition_id_t, comm_volume_t>;
  /// Partition adjacency list
  using part_adj_t = std::vector<part_edges_t>;

  /**
   *  This requires a network in which every vertex is assigned a partition.
   *  The number of partitions hints the allocation of internal vectors.
   */
  Partition_Info(const std::shared_ptr<const wcs::Network>& net,
                 const size_t num_partitions = 0ul,
                 const load_t min_load = load_t{0});

  void reset_num_partitions();
  void scan(bool verbose = false);
  void report() const;

 protected:
  /// Network with partition info
  std::shared_ptr<const wcs::Network> m_net;

  /// Number of partitions
  size_t m_num_parts;

  /// Minimum compute load that is assigned to the vertices with zero load
  load_t m_min_load;

  // The element v[t][p] of a vector v below stores the information of
  // type t vertices in partition p.

  /// Number of local vertices per type for each partition
  std::vector<p_sum_t> m_num_local;

  /// Sum of in-degrees of local vertices per vertex type per partition
  std::vector<p_sum_t> m_sum_indegree;

  /// Sum of out-degrees of local vertices per vertex type per partition
  std::vector<p_sum_t> m_sum_outdegree;

  /// Maximum in-degree of local vertices per vertex type per partition
  std::vector<p_max_t> m_max_indegree;

  /// Maximum out-degree of local vertices per vertex type per partition
  std::vector<p_max_t> m_max_outdegree;

  /// Maximum degree of local vertices per vertex type per partition
  std::vector<p_max_t> m_max_degree;

  /// Sum of rates of local reactions per vertex type per partition
  std::vector<p_load_t> m_load;

  /** Connectivity and outgoing communication volume of partition per type of
   *  neighboring vertex */
  std::vector<part_adj_t> m_comm_out;

  /** Connectivity and incoming communication volume of partition per type of
   *  neighboring vertex */
  std::vector<part_adj_t> m_comm_in;
};


/**@}*/
} // end of namespace wcs
#endif // __WCS_PARTITION_PARTITION_INFO_HPP__
