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

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

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
  using reaction_times_t = Sim_State_Change::reaction_times_t;

  using affected_reactions_t = Sim_State_Change::affected_reactions_t;
  using dependent_reactions_t = Sim_State_Change::dependent_reactions_t;


  SSA_NRM(const std::shared_ptr<wcs::Network>& net_ptr);
  SSA_NRM(SSA_NRM&& other) = default;
  SSA_NRM& operator=(SSA_NRM&& other) = default;
  ~SSA_NRM() override;

  void init(const sim_iter_t max_iter,
            const sim_time_t max_time,
            const unsigned rng_seed) override;

  /**
   * Determines the next reaction and the time when it occurs.
   * When successful, this function returns Sim_Method::Success. Otherwise,
   * it returns a failure code.
   */
  Sim_Method::result_t schedule(revent_t& evt);
  /**
   * Execute a reaction event.
   * Check the simulation termination condition at the beginning. If it is not
   * to be terminated yet, proceed and return true. Otherwise, stop immediately
   * and return false.
   */
  bool forward(const revent_t evt);
  /// Main loop of SSA
  std::pair<sim_iter_t, sim_time_t> run() override;

 #if defined(WCS_HAS_ROSS)
  void backward(revent_t& evt);
  void commit_des();
   
  /** Record as many states as the given number of iterations from the
   *  beginning of the digest list */
  void record_first_n(const sim_iter_t num) override;
 #endif // defined(WCS_HAS_ROSS)

  rng_t& rgen();

  /// Returns the earliest reaction
  priority_t choose_reaction() const;

  /// Check if priority the queue is empty
  bool is_empty() const;

 #if defined(_OPENMP) && defined(WCS_OMP_RUN_PARTITION)
  /// Used in parallel mode with partitioned network
  bool advance_time_and_iter(const sim_time_t t_new);
  void update_reactions(const sim_time_t t_fired,
                      #ifdef WCS_CACHE_DEPENDENT
                        const dependent_reactions_t& affected,
                      #else
                        const affected_reactions_t& affected,
                      #endif // WCS_CACHE_DEPENDENT
                        reaction_times_t& affected_rtimes);
 #endif // defined(_OPENMP) && defined(WCS_OMP_RUN_PARTITION)

  void update_reactions(const priority_t& fired,
                      #ifdef WCS_CACHE_DEPENDENT
                        const dependent_reactions_t& affected,
                      #else
                        const affected_reactions_t& affected,
                      #endif // WCS_CACHE_DEPENDENT
                        reaction_times_t& affected_rtimes);

protected:
  void build_heap();
  sim_time_t get_reaction_time();
  wcs::sim_time_t recompute_reaction_time(const v_desc_t& vd);
  wcs::sim_time_t adjust_reaction_time(const v_desc_t& vd, wcs::sim_time_t rt);
  void revert_reaction_updates(const reaction_times_t& affected);

  void save_rgen_state(Sim_State_Change& digest) const;
  void load_rgen_state(const Sim_State_Change& digest);

protected:
  /** In-heap index table maintains where in the heap each item can be found.
      The position in the heap is identified by an index of type idx_t. */
  using heap_idx_t = int; // the type used in iheap
  using in_heap_index_table_t = std::unordered_map<v_desc_t, heap_idx_t>;

  in_heap_index_table_t m_idx_table;
  priority_queue_t m_heap;
  rng_t m_rgen;

 #if defined(WCS_HAS_ROSS)
  using digest_list_t = std::list<Sim_State_Change>;
  digest_list_t m_digests;
 #endif // defined(WCS_HAS_ROSS)
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_SIM_METHODS_SSA_NRM_HPP__
