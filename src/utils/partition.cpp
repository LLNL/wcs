/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include <unordered_map>
#include <type_traits>
#include <vector>

#include <utils/partition.hpp>
#include <utils/detect_methods.hpp>
#include <utils/exception.hpp>

namespace wcs {
/** \addtogroup wcs_utils
 *  *  @{ */

void Partition::make_intermediate_graph(const wcs::Network& net)
{
  using G_org = wcs::Network::graph_t; // The type of the original graph
  const G_org& g_org = net.graph(); // The original graph

  // Make sure if the original graph is bidirectional such that we can easily
  // find the incoming edges as well as the outgoing ones of a reaction vertex
  using directed_org = typename boost::graph_traits<G_org>::directed_category;
  constexpr bool is_org_bidirectional
    = std::is_same<directed_org, boost::bidirectional_tag>::value;

  if constexpr (!is_org_bidirectional) {
    WCS_THROW("The original graph is not bidirectional!");
    return;
  }

  // Edge property type of the new graph
  using ep_new_t = typename boost::edge_bundle_type<intm_graph_t>::type;

  // Make sure the new graph supports the weighted edge
  if constexpr (!has_get_weight<ep_new_t>::value ||
                !has_set_weight<ep_new_t>::value)
  {
    WCS_THROW("The new graph does not have get/set_weight()!");
    return;
  }

  // Reserve space for vertices in the new graph
  if constexpr (has_reserve_for_vertex_list<intm_graph_t>::value) {
    m_g_intm.m_vertices.reserve(net.get_num_species());
  }

  // The type of the vertex descriptor of the original graph
  using vd_org_t = wcs::Network::v_desc_t;

  // The type of the vertex descriptor of the new graph
  using vd_new_t = typename boost::graph_traits<intm_graph_t>::vertex_descriptor;

  // Map from the vertex descriptor of the original graph to that of the new
  // copy added to the new graph
  typename std::unordered_map<vd_org_t, vd_new_t> vd_org2new;

  // List of the descriptors of reaction vertices
  std::vector<vd_org_t> reactions;
  reactions.reserve(net.get_num_reactions()); // Reserve the space

  // Add vertices to the new graph
  typename boost::graph_traits<G_org>::vertex_iterator vi, vi_end;

  for (boost::tie(vi, vi_end) = boost::vertices(g_org); vi != vi_end; ++vi) {
    using vp_org_t = wcs::Network::v_prop_t;
    const vp_org_t& vp = g_org[*vi];
    const auto vt = static_cast<vp_org_t::vertex_type>(vp.get_typeid());

    if (vt == wcs::Vertex::_reaction_) { // Save the reaction vertex and move on
      reactions.push_back(*vi);
      continue;
    }

    // Make a copy of the non-reaction vertex and add it to the new graph
    vd_new_t vd = boost::add_vertex(wcs::Vertex{vp}, m_g_intm);

    // Add a mapping from the original vertex descriptor to the new one
    vd_org2new[*vi] = vd;
  }

  if (reactions.empty()) {
    using std::operator<<;
    std::cout << "Either there is no reaction in the network or "
              << "it has not been initialized by calling Network::init()"
              << std::endl;
    return;
  }

  // Add edges to the new graph
  for (const auto& reaction : reactions) {
    // Property of this reaction vertex
    auto & rp = g_org[reaction].property< Reaction<vd_org_t> >();
    // The rate of this reaction
    const auto rate = rp.get_rate();

    // Outputs of the reaction, such as product species 
    std::vector<vd_new_t> outputs;
    // Inputs of the reaction and its rate formula, such as the reactants
    std::vector<vd_new_t> inputs;

    // Collect the vertex descriptors in the new graph corresponding to
    // the output vertices of this reaction
    for(const auto ei_out :
        boost::make_iterator_range(boost::out_edges(reaction, g_org))) {
      outputs.emplace_back(vd_org2new.at(boost::target(ei_out, g_org)));
    }

    // Collect the vertex descriptors in the new graph corresponding to
    // the input vertices of this reaction
    for(const auto ei_in :
        boost::make_iterator_range(boost::in_edges(reaction, g_org))) {
      inputs.emplace_back(vd_org2new.at(boost::source(ei_in, g_org)));
    }

    // The number of edges. If no input is reused as an output,
    // it is (num of inputs) * (num of outputs)
    size_t num_edges = 0ul;

    for(const auto& v_in : inputs) {
      bool reappear = false;
      for (const auto& v_out : outputs) {
        reappear |= (v_in == v_out);
      }
      if (reappear) { // Don't consider an edge from one to itself
        num_edges += outputs.size() - 1ul;
      } else {
        num_edges += outputs.size();
      }
    }


    // Create edge between inputs and outputs of the reaction
    for(const auto& v_in : inputs) {
      for (const auto& v_out : outputs) {
        if (v_in == v_out) {
          continue;
        }
        auto ret = boost::add_edge(v_in, v_out, m_g_intm);
        if (!ret.second) {
          WCS_THROW("Failed to add an edge from an input to an output.");
          return;
        }
        auto& e_prop = m_g_intm[ret.first];
        e_prop.set_weight(rate/static_cast<wcs::reaction_rate_t>(num_edges));
        // An edge is never cut in the vetex-cut partitioning. Therefore, the
        // weight of an edge does not encode the penalty for splitting it
        // over two partitions. Instead, it represent the compute load to
        // balance.
      }
    }

    // Create an edge between the inputs to make sure that all the inputs
    // and the outputs involved in the reaction stay in the same partition
    // after vertex-cut partitioning.
    for(size_t i = 1ul; i < inputs.size(); ++i) {
      auto ret = boost::add_edge(inputs[0], inputs[i], m_g_intm);
      if (!ret.second) {
        WCS_THROW("Failed to add an edge from an input to another input.");
        return;
      }
      auto& e_prop = m_g_intm[ret.first];
      e_prop.set_weight(static_cast<wcs::reaction_rate_t>(0));
    }
  }
}

const Partition::intm_graph_t& Partition::intm_graph() const
{
  return m_g_intm;
}

Partition::intm_graph_t& Partition::intm_graph()
{
  return m_g_intm;
}

/**@}*/
} // end of namespace wcs
