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
#include "reaction_network/network.hpp"

namespace wcs {
/** \addtogroup wcs_sim_methods
 *  @{ */

struct Sim_State_Change {
  using v_desc_t = wcs::Network::v_desc_t;

  /** Type for keeping track of species updates to facilitate undoing
   *  reaction processing.  */
  using update_t = std::pair<v_desc_t, stoic_t>;
  using update_list_t = std::vector<update_t>;
  /** Type for the list of reactions that share any of the species with the
   *  firing reaction */
  using affected_reactions_t = std::set<v_desc_t>;
  using reaction_times_t = std::vector<std::pair<v_desc_t, sim_time_t> >;

  sim_time_t m_sim_time; ///< Current simulation time
  v_desc_t m_reaction_fired;

  /// species updated as a result of firing the reaction
  update_list_t m_species_updated;
  /**
   *  Any other reaction that takes a reactant species being updated as a
   *  result of processing the reaction.
   */
  affected_reactions_t m_reactions_affected;

  /**
   *  The currently expected times of the affected reactions to occur.
   *  This is used by the next reaction method.
   */
  reaction_times_t m_reaction_times;

  /// Serialized RNG states
  std::vector<char> m_rng_state;

  void clear() {
    m_sim_time = wcs::Network::get_etime_ulimit();
    m_species_updated.clear();
    m_reactions_affected.clear();
    m_reaction_times.clear();
    m_rng_state.clear();
  }
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_SIM_METHODS_SIM_STATE_CHANGE_HPP__
