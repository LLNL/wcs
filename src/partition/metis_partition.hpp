/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef  __WCS_PARTITION_METIS_PARTITION__
#define  __WCS_PARTITION_METIS_PARTITION__

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

// Only enable when METIS is available
#if defined(WCS_HAS_METIS)
#include <metis.h>
#include <vector>
#include <map>
#include <iostream>
#include <unordered_map>
#include <type_traits> // is_same
#include "reaction_network/network.hpp"
#include "partition/metis_params.hpp"

namespace wcs {
/** \addtogroup wcs_partition
 *  @{ */

class Metis_Partition {
 public:
  using graph_t = wcs::Network::graph_t;
  using directed_category = typename boost::graph_traits<graph_t>::directed_category;
  using v_desc_t = wcs::Network::v_desc_t;

  Metis_Partition(const Metis_Params& mp);

  /** Return the number of undirected edges that can go into the header of
    * a Metis input file */
  size_t get_num_edges() const;
  /// Return the number of vertices
  size_t get_num_vertices() const;

  /// Populate input graph data structure for Metis
  void prepare();
  /// Check Metis return code
  static bool check_run(const int ret, const bool verbose = false);
  /// Run partition
  bool run(std::vector<idx_t>& parts, idx_t& objval);

  /// Print partitioning parameters
  void print_params() const;

  /// Print Metis input graph (using vertext index)
  static void print_metis_inputs(const std::vector<idx_t>& xadj,
                                 const std::vector<idx_t>& adjcny,
                                 const std::vector<idx_t>& vwgt,
                                 const std::vector<idx_t>& vsize,
                                 std::ostream& os = std::cout);

  /// Print the Metis input graph (using vertext index) of this object
  void print_metis_graph(std::ostream& os = std::cout) const;

  /// Print input adjacency list of the graph (using vertex name)
  void print_adjacency(std::ostream& os = std::cout) const;


  const std::unordered_map<v_desc_t, idx_t>& get_map_from_desc_to_idx() const;
  const std::map<idx_t, v_desc_t>& get_map_from_idx_to_desc() const;

 protected:
  using v_iter_t = wcs::Network::v_iter_t;
  using v_prop_t = wcs::Network::v_prop_t;
  using r_prop_t = wcs::Network::r_prop_t;
  using reaction_rate_t = wcs::reaction_rate_t;

  /**
   * Build a map from the BGL vertex descriptor to the sequential index of
   * idx_t type used in Metis, and find out the total number of edges.
   */
  void build_map_from_desc_to_idx();
  /// Populate the adjacency list in compressed storage format for Metis
  void populate_adjacny_list();
  /// Populate the list of vertex weights and the list of vertex sizes
  void populate_vertex_info();

 protected:
  static constexpr bool is_bidirectional
    = std::is_same<directed_category, boost::bidirectional_tag>::value;

  /// Metis configuration parameters
  Metis_Params m_p;

  // idx_t is the type defined in metis

  /// Map from the BGL vertex descriptor to the vertex index used in Metis
  std::unordered_map<v_desc_t, idx_t> m_vd2idx;
  /// Map from the vertex index used in Metis to the BGL vertex descriptor
  std::map<idx_t, v_desc_t> m_idx2vd;

  /// Indices to the adjacency per vertex in the CSR format
  std::vector<idx_t> m_xadj;
  /// Adjacency list in CSR format
  std::vector<idx_t> m_adjncy;

  /// Vertex weights (the list of load balancing constraint sets for each vertex)
  std::vector<idx_t> m_vwgt;
  /// Vertex sizes used in computing communication volume
  std::vector<idx_t> m_vsize;

  /// Number of undirected edges recognized by Metis partitioning
  size_t m_num_edges;
  /// Number of vertices in the graph
  size_t m_num_vertices;
};

/**@}*/
} // end of namespace wcs
#endif // defined(WCS_HAS_CONFIG)
#endif //  __WCS_PARTITION_METIS_PARTITION__
