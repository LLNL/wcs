/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __WCS_SIM_METHODS_SIM_STATE_CHANGE_HPP__
#define __WCS_SIM_METHODS_SIM_STATE_CHANGE_HPP__
#include <vector>
#include "wcs_types.hpp"
#include "sim_methods/update.hpp"
#include "reaction_network/network.hpp"

//#define ENABLE_SPECIES_UPDATE_TRACKING

namespace wcs {
/** \addtogroup wcs_sim_methods
 *  @{ */

struct Sim_State_Change {
  using v_desc_t = wcs::Network::v_desc_t;

  /** Type for keeping track of species updates to facilitate undoing
   *  reaction processing.  */
  /** Type for the list of reactions that share any of the species with the
   *  firing reaction */
  using affected_reactions_t = std::set<v_desc_t>;
  /** To carry the same information as affected_reactions_t does. However,
   *  by using the randomly accessible container, this makes it easier to
   *  construct a parallel loop over the container.
   */
  using dependent_reactions_t = std::vector<v_desc_t>;
  using reaction_times_t = std::vector<std::pair<v_desc_t, sim_time_t> >;
  using revent_t = std::pair<sim_time_t, v_desc_t>;

  sim_time_t m_sim_time; ///< Current simulation time
  v_desc_t m_reaction_fired;

 #ifdef ENABLE_SPECIES_UPDATE_TRACKING
  /// species counts updated as a result of firing the reaction
  cnt_updates_t m_species_updated;
 #endif // ENABLE_SPECIES_UPDATE_TRACKING

  /// species concentrations updated as a result of firing the reaction
  //conc_updates_t m_species_updated;

  /**
   *  Any other reaction that takes a reactant species being updated as a
   *  result of processing the reaction.
   */
  affected_reactions_t m_reactions_affected;
  dependent_reactions_t m_dependent_reactions;

  /**
   *  The currently expected times of the affected reactions to occur.
   *  This is used by the next reaction method.
   */
  reaction_times_t m_reaction_times;

  /// Serialized RNG states
  std::vector<char> m_rng_state;

  Sim_State_Change(const revent_t& firing)
  : m_sim_time(firing.first), m_reaction_fired(firing.second) {}

  Sim_State_Change() = default;
  Sim_State_Change(const Sim_State_Change& other) = default;
  Sim_State_Change(Sim_State_Change&& other) = default;
  Sim_State_Change& operator=(const Sim_State_Change& other) = default;
  Sim_State_Change& operator=(Sim_State_Change&& other) = default;

  void clear() {
    m_sim_time = wcs::Network::get_etime_ulimit();
   #ifdef ENABLE_SPECIES_UPDATE_TRACKING
    m_species_updated.clear();
   #endif // ENABLE_SPECIES_UPDATE_TRACKING
    m_reactions_affected.clear();
    m_reaction_times.clear();
    m_rng_state.clear();
  }
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_SIM_METHODS_SIM_STATE_CHANGE_HPP__
