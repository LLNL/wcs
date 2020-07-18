/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __WCS_SIM_METHODS_SSA_NRM_HPP__
#define __WCS_SIM_METHODS_SSA_NRM_HPP__
#include <cmath>
#include <limits>
#include "sim_methods/sim_method.hpp"

namespace wcs {
/** \addtogroup wcs_sim_methods
 *  @{ */

class SSA_NRM : public Sim_Method {
public:
  using rng_t = wcs::RNGen<std::uniform_int_distribution, unsigned>;
  using v_desc_t = Sim_Method::v_desc_t;
  using priority_t = std::pair<wcs::sim_time_t, v_desc_t>;
  /// Type of heap structure
  using priority_queue_t = std::vector<priority_t>;
  /** Type for the list of reactions that share any of the species with the
   *  firing reaction */
  using affected_reactions_t = std::set<v_desc_t>;
  /// Type of the pair of BGL vertex descriptor for reaction and the its time
  using reaction_times_t = std::vector<std::pair<v_desc_t, wcs::sim_time_t> >;

  /** Type for keeping track of species updates to facilitate undoing
   *  reaction processing.  */
  using update_t = std::pair<v_desc_t, stoic_t>;
  using update_list_t = std::vector<update_t>;

  SSA_NRM();
  ~SSA_NRM() override;
  void init(std::shared_ptr<wcs::Network>& net_ptr,
            const sim_iter_t max_iter,
            const sim_time_t max_time,
            const unsigned rng_seed) override;

  std::pair<sim_iter_t, sim_time_t> run() override;

  static bool later(const priority_t& v1, const priority_t& v2);
  rng_t& rgen();

protected:
  void build_heap();
  priority_t& choose_reaction();
  sim_time_t get_reaction_time(const priority_t& p);
  bool fire_reaction(const v_desc_t vd_firing,
                     update_list_t& updating_species,
                     affected_reactions_t& affected_reactions);
  void reset_reaction_time(const v_desc_t& vd, wcs::sim_time_t& rt);
  void adjust_reaction_time(const v_desc_t& vd, wcs::sim_time_t& rt);
  void update_reactions(priority_t& firing,
                       const affected_reactions_t& affected,
                       reaction_times_t& affected_rtimes);
  void revert_reaction_updates(const sim_time_t dt,
                       const reaction_times_t& affected);
  void undo_species_updates(const update_list_t& updates) const;
  bool undo_reaction(const v_desc_t vd_undo,
                     update_list_t& reverting_species,
                     affected_reactions_t& affected_reactions);

protected:
  priority_queue_t m_heap;
  rng_t m_rgen;
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_SIM_METHODS_SSA_NRM_HPP__
