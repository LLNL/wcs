/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __WCS_REACTION_NETWORK_NETWORK_HPP__
#define __WCS_REACTION_NETWORK_NETWORK_HPP__
// Whole Cell Model Simulator
/** @file
 * \ingroup wcs_reaction_network
 * \brief reaction graph representation.
 */
/** \ingroup wcs_reaction_network
 * \class wcs::Network
 * \brief reaction graph representation as a bipartite graph.
 *
 * Represent a reaction network as a bipartite graph using BGL.
 * A reaction graph consists of a set of species nodes, a set of reaction nodes,
 * and the directed edges between species and reactions.
 *
 * \author Jae-Seung Yeom <yeom2@llnl.gov>
 * \date 2019
 */

#include <string>
#include <unordered_map>
#include "bgl.hpp"
#include "reaction_network/species.hpp"
#include "reaction_network/reaction.hpp"
#include "reaction_network/vertex.hpp"
#include "reaction_network/edge.hpp"

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  *  @{ */

class Network {
 public:
    /// The type of the BGL graph to represent reaction networks
  using graph_t  = boost::adjacency_list<
                   wcs_out_edge_list_t,
                   wcs_vertex_list_t,
                   boost::bidirectionalS,
                   wcs::Vertex, // vertex property bundle
                   wcs::Edge,   // edge property bundle
                   boost::no_property,
                   boost::vecS>;

  /// The type of the vertex property bundle
  using v_prop_t = boost::vertex_bundle_type<graph_t>::type;
  /// The type of BGL vertex descriptor for graph_t
  using v_desc_t = boost::graph_traits<graph_t>::vertex_descriptor;
  /// The type of BGL vertex iterator for graph_t
  using v_iter_t = boost::graph_traits<graph_t>::vertex_iterator;
  /// The type of the edge property bundle
  using e_prop_t = boost::edge_bundle_type<graph_t>::type;
  /// The type of BGL edge descriptor for graph_t
  using e_desc_t = boost::graph_traits<graph_t>::edge_descriptor;
  /// The type of BGL edge iterator for graph_t
  using e_iter_t = boost::graph_traits<graph_t>::edge_iterator;

  using vertex_type = v_prop_t::vertex_type;
  using v_label_t = std::string;

  // Type of the species name
  using s_label_t = v_label_t;
  // reaction driver type, std::pair<v_desc_t, stoic_t>
  using rdriver_t = Reaction<v_desc_t>::rdriver_t;
  // Type of the map of species involve in a reaction
  using s_involved_t = std::map<s_label_t, rdriver_t>;
  // Reaction descriptor type
  using r_desc_t = std::pair<v_desc_t, s_involved_t>;

  using rand_access
    = typename boost::detail::is_random_access<typename graph_t::vertex_list_selector>::type;

  using reaction_list_t = std::vector<v_desc_t>;

  /** The type of the list of species. This is chosen for the memory efficiency
      than for the lookup performance. */
  using species_list_t  = std::vector<v_desc_t>;

 public:
  /// Load an input Graph ML file
  void load(const std::string graphml_filename);
  void init();
  reaction_rate_t set_reaction_rate(const v_desc_t r) const;
  reaction_rate_t get_reaction_rate(const v_desc_t r) const;

  size_t get_num_vertices() const;
  size_t get_num_species() const;
  size_t get_num_reactions() const;
  size_t get_num_vertices_of_type(const vertex_type vt) const;

  /// Allow read-only access to the internal BGL graph
  const graph_t& graph() const;
  /// Allow read-only access to the internal reaction list
  const reaction_list_t& reaction_list() const;
  /// Allow read-only access to the internal species list
  const species_list_t& species_list() const;
  /// Find the species by the label and return the BGL vertex descriptor
  v_desc_t find_species(const std::string& label);
  /// Set the largest delay period for an active reaction to fire
  static void set_etime_ulimit(const sim_time_t t);
  /// Return the largest delay period for an active reaction to fire
  static sim_time_t get_etime_ulimit();

  bool check_reaction(const v_desc_t r) const;

  std::string show_species_labels(const std::string title="Species: ") const;
  std::string show_reaction_labels(const std::string title="Reaction:") const;
  std::string show_species_counts() const;
  std::string show_reaction_rates() const;

 protected:
  /// Sort the species list by the label (in lexicogrphical order)
  void sort_species();
  void loadGraphML(const std::string graphml_filename);
  void loadSBML(const std::string sbml_filename);

 protected:
  /// The BGL graph to represent a reaction network
  graph_t m_graph;

  /// List of the BGL descriptors of reaction type vertices
  reaction_list_t m_reactions;

  /// List of the BGL descriptors of species type vertices
  species_list_t m_species;

  /**
   * The upper limit of the delay period for an active reaction to fire beyond
   * which we consider the reaction inactive/disabled. This is by default set to
   * the largest value of the sim_time_t type.
   */
  static sim_time_t m_etime_ulimit;
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_REACTION_NETWORK_NETWORK_HPP__
