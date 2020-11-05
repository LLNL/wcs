/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __WCS_SIM_METHODS_SSA_DIRECT_HPP__
#define __WCS_SIM_METHODS_SSA_DIRECT_HPP__

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

/**
 *  This direct SSA method implementation takes advantage of the dependency
 *  graph to precisely identify the propensities to update and performs the
 *  binary search over a propensity list sorted in the ascending order.
 *  This differs from the Cao's optimized direct method where the propensity
 *  list is ordered in the descending order to reduce the average linear
 *  search length.
 *  This reverses the order to minimize the potential numerical error caused
 *  by adding a small propensity to a large cumulative propensity when
 *  searching a value from the largest propensity to the smallest propensity.
 */
class SSA_Direct : public Sim_Method {
public:
  using rng_t = wcs::RNGen<std::uniform_real_distribution, double>;
  using v_desc_t = Sim_Method::v_desc_t;
  using priority_t = std::pair<reaction_rate_t, v_desc_t>;
  using propensisty_list_t = std::vector<priority_t>;


  SSA_Direct(const std::shared_ptr<wcs::Network>& net_ptr);
  SSA_Direct(SSA_Direct&& other) = default;
  SSA_Direct& operator=(SSA_Direct&& other) = default;
  ~SSA_Direct() override;

  /// Initialize propensity list
  void init(const unsigned max_iter,
            const double max_time,
            const unsigned rng_seed) override;

  /// Main loop of SSA
  std::pair<unsigned, sim_time_t> run() override;
  Sim_Method::result_t forward(Sim_State_Change& digest);

 #if defined(WCS_HAS_ROSS)
  void backward(Sim_State_Change& digest);

  /** Record as many states as the given number of iterations from the
   *  beginning of the digest list */
  void record_first_n(const sim_iter_t num) override;
 #endif // defined(WCS_HAS_ROSS)

  static bool less(const priority_t& v1, const priority_t& v2);

  rng_t& rgen_e();
  rng_t& rgen_t();

protected:
  void build_propensity_list();
  priority_t& choose_reaction();
  sim_time_t get_reaction_time();
  void update_reactions(priority_t& fired,
                        const Sim_Method::affected_reactions_t& affected,
                        const bool check_reaction);

  void save_rgen_state(Sim_State_Change& digest);
  void load_rgen_state(const Sim_State_Change& digest);
  Sim_Method::result_t schedule();

protected:
  /// Cumulative propensity of reactions events
  propensisty_list_t m_propensity;
  rng_t m_rgen_evt; ///< RNG for events
  rng_t m_rgen_tm; ///< RNG for event times
  /// map from vertex descriptor to propensity
  std::unordered_map<v_desc_t, size_t> m_pindices;

 #if defined(WCS_HAS_ROSS)
  using digest_list_t = std::list<Sim_State_Change>;
  digest_list_t m_digests;
 #endif // defined(WCS_HAS_ROSS)
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_SIM_METHODS_SSA_DIRECT_HPP__
