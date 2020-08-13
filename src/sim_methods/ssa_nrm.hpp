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
#include <unordered_map>
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
  /// Type of the pair of BGL vertex descriptor for reaction and the its time
  using reaction_times_t = std::vector<std::pair<v_desc_t, wcs::sim_time_t> >;


  SSA_NRM();
  ~SSA_NRM() override;
  void init(std::shared_ptr<wcs::Network>& net_ptr,
            const sim_iter_t max_iter,
            const sim_time_t max_time,
            const unsigned rng_seed) override;

  std::pair<sim_iter_t, sim_time_t> run() override;

  rng_t& rgen();

protected:
  void build_heap();
  priority_t choose_reaction();
  sim_time_t get_reaction_time(const priority_t& p);
  wcs::sim_time_t recompute_reaction_time(const v_desc_t& vd);
  wcs::sim_time_t adjust_reaction_time(const v_desc_t& vd, wcs::sim_time_t rt);
  void update_reactions(const priority_t& fired,
                        const Sim_Method::affected_reactions_t& affected,
                        reaction_times_t& affected_rtimes);
  void revert_reaction_updates(const sim_time_t dt,
                               const reaction_times_t& affected);

protected:
  /** In-heap index table maintains where in the heap each item can be found.
      The position in the heap is identified by an index of type idx_t. */
  using heap_idx_t = int; // the type used in iheap
  using in_heap_index_table_t = std::unordered_map<v_desc_t, heap_idx_t>;

  in_heap_index_table_t m_idx_table;
  priority_queue_t m_heap;
  rng_t m_rgen;
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_SIM_METHODS_SSA_NRM_HPP__
