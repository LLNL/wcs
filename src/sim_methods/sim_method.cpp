/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include "sim_methods/sim_method.hpp"
#include <limits>

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  @{ */

Sim_Method::Sim_Method()
: m_net_ptr(nullptr),
  m_max_iter(static_cast<sim_iter_t>(0u)),
  m_max_time(static_cast<sim_time_t>(0)),
  m_cur_iter(static_cast<sim_iter_t>(0u)),
  m_sim_time(static_cast<sim_time_t>(0)),
  m_enable_tracing(false),
  m_enable_sampling(false),
  m_sample_iter_interval(static_cast<sim_iter_t>(0u)),
  m_sample_time_interval(static_cast<sim_time_t>(0)),
  m_next_sample_iter(static_cast<sim_iter_t>(0u)),
  m_next_sample_time(static_cast<sim_time_t>(0))
{
  using directed_category
    = typename boost::graph_traits<wcs::Network::graph_t>::directed_category;

  constexpr bool is_bidirectional
    = std::is_same<directed_category, boost::bidirectional_tag>::value;

  if constexpr (!is_bidirectional) {
    WCS_THROW("Cannot get species population without in-edges.");
  }
}

Sim_Method::~Sim_Method() {}

void Sim_Method::set_tracing()
{
  m_enable_tracing = true;
  m_enable_sampling = false;
}

void Sim_Method::unset_tracing()
{
  m_enable_tracing = false;
}

void Sim_Method::set_sampling(const sim_time_t time_interval)
{
  m_enable_tracing = false;
  m_enable_sampling = true;
  m_sample_time_interval = time_interval;
  m_next_sample_time = m_sim_time + time_interval;
  m_next_sample_iter = std::numeric_limits<sim_iter_t>::max();
}

void Sim_Method::set_sampling(const sim_iter_t iter_interval)
{
  m_enable_tracing = false;
  m_enable_sampling = true;
  m_sample_iter_interval = iter_interval;
  m_next_sample_iter = m_cur_iter + iter_interval;
  m_next_sample_time = std::numeric_limits<sim_time_t>::infinity();
}

void Sim_Method::unset_sampling()
{
  m_enable_sampling = false;
}

void Sim_Method::record_initial_state(const std::shared_ptr<wcs::Network>& net_ptr)
{ // record initial state of the network
  if (m_enable_tracing) {
    m_trace.record_initial_condition(m_net_ptr);
  } else if (m_enable_sampling) {
    m_samples.record_initial_condition(m_net_ptr);
  }
}

void Sim_Method::record_final_state(const v_desc_t rv)
{
  if (m_enable_tracing) {
    m_trace.record_reaction(m_sim_time, rv);
  } else if (m_enable_sampling) {
    m_samples.record_reaction(rv);
    m_samples.take_sample(m_sim_time);
  }
}

bool Sim_Method::check_to_record(const v_desc_t rv)
{
  if (m_enable_tracing) {
    m_trace.record_reaction(m_sim_time, rv);
    return true;
  } else if (m_enable_sampling) {
    m_samples.record_reaction(rv);
    if (m_cur_iter >= m_next_sample_iter) {
      m_next_sample_iter += m_sample_iter_interval;
      m_samples.take_sample(m_sim_time);
      return true;
    } else if (m_sim_time >= m_next_sample_time) {
      m_next_sample_time += m_sample_time_interval;
      m_samples.take_sample(m_sim_time);
      return true;
    } else {
      return false;
    }
  }
  return false;
}

void Sim_Method::pop_trace()
{
  if (m_enable_tracing) {
    m_trace.pop_back();
  }
}

bool Sim_Method::check_to_record()
{
  if (m_enable_tracing) {
    return true;
  } else if (m_enable_sampling) {
    if (m_cur_iter >= m_next_sample_iter) {
      m_next_sample_iter += m_sample_iter_interval;
      return true;
    } else if (m_sim_time >= m_next_sample_time) {
      m_next_sample_time += m_sample_time_interval;
      return true;
    } else {
      return false;
    }
  }
  return false;
}

Sim_Method::trace_t& Sim_Method::trace() {
  return m_trace;
}

Sim_Method::samples_t& Sim_Method::samples() {
  return m_samples;
}


/**
 * Execute the chosen reaction.
 * In other words, update the species population involved in the reaction.
 * In addition, record how the species are updated, and which other reactions
 * are affected as a result.
 * The former can be used to undo the reaction if needed. The latter is used
 * to update the propensity of the reactions affected by the changes in species
 * counts.
 */
bool Sim_Method::fire_reaction(Sim_State_Change& digest)
{
  using s_prop_t = wcs::Species;

  // Reactions do not change the connectivity, but only change the property
  // data, which are allocated outside of BGL graph, but only linked to it.
  const wcs::Network::graph_t& g = m_net_ptr->graph();
  // The vertex descriptor of the reaction to undo
  const auto& rd_firing = digest.m_reaction_fired;
  auto& updating_species = digest.m_species_updated;
  auto& reactions_affected = digest.m_reactions_affected;

  updating_species.clear();
  reactions_affected.clear();

  // reactant species
  for (const auto ei_in :
       boost::make_iterator_range(boost::in_edges(rd_firing, g)))
  {
    const auto vd_updating = boost::source(ei_in, g);
    const auto& sv_updating = g[vd_updating];
    if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
      // in case that there are other type of vertices than species or reaction
      if (sv_updating.get_type() != wcs::Vertex::_species_) continue;
    }

    auto& sp_updating = sv_updating.property<s_prop_t>();
    const auto stoichio = g[ei_in].get_stoichiometry_ratio();
    if (stoichio == static_cast<stoic_t>(0)) {
      continue;
    }
  #ifdef NDEBUG
    sp_updating.dec_count(stoichio);
  #else
    // This really should not happen because whether the reaction is feasible is
    // checked before computing reaction time or propensity.
    if (!sp_updating.dec_count(stoichio)) { // State update
      std::string err = "Not enough reactants of " + sv_updating.get_label()
                      + "[" + std::to_string(sp_updating.get_count())
                      + "] for reaction " + g[rd_firing].get_label();
      WCS_THROW(err);
      return false;
    }
  #endif
    updating_species.emplace_back(std::make_pair(vd_updating, -stoichio));

    for (const auto vi_affected :
         boost::make_iterator_range(boost::out_edges(vd_updating, g)))
    {
      const auto rd_affected = boost::target(vi_affected, g);
      if (rd_affected == rd_firing) continue;
      reactions_affected.insert(rd_affected);
    }
  }

  // product species
  for (const auto ei_out :
       boost::make_iterator_range(boost::out_edges(rd_firing, g)))
  {
    const auto vd_updating = boost::target(ei_out, g);
    const auto& sv_updating = g[vd_updating];
    if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
      // in case that there are other type of vertices than species or reaction
      if (sv_updating.get_type() != wcs::Vertex::_species_) continue;
    }

    auto& sp_updating = sv_updating.property<s_prop_t>();
    const auto stoichio = g[ei_out].get_stoichiometry_ratio();
    if (stoichio == static_cast<stoic_t>(0)) {
      continue;
    }
  #ifdef NDEBUG
    sp_updating.inc_count(stoichio);
  #else
    if (!sp_updating.inc_count(stoichio)) { // State update
      std::string err = "Can not produce more of " + sv_updating.get_label()
                      + "[" + std::to_string(sp_updating.get_count())
                      + "] by reaction " + g[rd_firing].get_label();
      WCS_THROW(err);
      return false;
    }
  #endif
    updating_species.emplace_back(std::make_pair(vd_updating, stoichio));

    for (const auto vi_affected :
         boost::make_iterator_range(boost::out_edges(vd_updating, g)))
    {
      const auto rd_affected = boost::target(vi_affected, g);
      if (rd_affected == rd_firing) continue;
      reactions_affected.insert(rd_affected);
    }
  }

  return true;
}


/**
 * Undo the species updates applied during incomplete reaction processing.
 * This relies on the list of updates made to species, and revert them.
 */
void Sim_Method::undo_species_updates(
  const Sim_Method::update_list_t& updates) const
{
  bool ok = true;
  const wcs::Network::graph_t& g = m_net_ptr->graph();

  for (const auto& u: updates) {
    const auto& sv_undo = g[u.first];
    using s_prop_t = wcs::Species;
    auto& sp_undo = sv_undo.property<s_prop_t>();
    if (u.second > static_cast<stoic_t>(0)) {
      ok &= sp_undo.dec_count(u.second);
    } else {
      ok &= sp_undo.inc_count(u.second);
    }
    if (!ok) {
      WCS_THROW("Failed to reverse the species updates");
    }
  }
}

/**
 * Undo the species states updated by a reaction identified by the given
 * vertex descriptor.
 * It is similar to `fire_reaction()` except that it changes the species count
 * in an opposite way.
 */
bool Sim_Method::undo_reaction(const Sim_Method::v_desc_t& rd_undo) const
{
  using s_prop_t = wcs::Species;

  // Reactions do not change the connectivity, but only change the property
  // data, which are allocated outside of BGL graph, but only linked to it.
  const wcs::Network::graph_t& g = m_net_ptr->graph();

  // reactant species
  for (const auto ei_in :
       boost::make_iterator_range(boost::in_edges(rd_undo, g)))
  {
    const auto vd_reverting = boost::source(ei_in, g);
    const auto& sv_reverting = g[vd_reverting];
    if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
      // in case that there are other type of vertices than species or reaction
      if (sv_reverting.get_type() != wcs::Vertex::_species_) continue;
    }

    auto& sp_reverting = sv_reverting.property<s_prop_t>();
    const auto stoichio = g[ei_in].get_stoichiometry_ratio();
    if (stoichio == static_cast<stoic_t>(0)) {
      continue;
    }
  #ifdef NDEBUG
    sp_reverting.inc_count(stoichio);
  #else
    if (!sp_reverting.inc_count(stoichio)) { // State update
      std::string err = "Unable to undo the decrement of reactant "
                      + sv_reverting.get_label()
                      + "[" + std::to_string(sp_reverting.get_count())
                      + "] for reaction " + g[rd_undo].get_label();
      WCS_THROW(err);
      return false;
    }
  #endif
  }

  // product species
  for (const auto ei_out :
       boost::make_iterator_range(boost::out_edges(rd_undo, g)))
  {
    const auto vd_reverting = boost::target(ei_out, g);
    const auto& sv_reverting = g[vd_reverting];
    if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
      // in case that there are other type of vertices than species or reaction
      if (sv_reverting.get_type() != wcs::Vertex::_species_) continue;
    }

    auto& sp_reverting = sv_reverting.property<s_prop_t>();
    const auto stoichio = g[ei_out].get_stoichiometry_ratio();
    if (stoichio == static_cast<stoic_t>(0)) {
      continue;
    }
  #ifdef NDEBUG
    sp_reverting.dec_count(stoichio);
  #else
    if (!sp_reverting.dec_count(stoichio)) { // State update
      std::string err = "Unable to undo the production of "
                      + sv_reverting.get_label()
                      + "[" + std::to_string(sp_reverting.get_count())
                      + "] by reaction " + g[rd_undo].get_label();
      WCS_THROW(err);
      return false;
    }
  #endif
  }

  return true;
}
/**@}*/
} // end of namespace wcs
