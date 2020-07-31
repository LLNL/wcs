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
#include <cmath>
#include <limits>
#include <unordered_map>
#include "sim_methods/sim_method.hpp"

namespace wcs {
/** \addtogroup wcs_sim_methods
 *  @{ */

class SSA_Direct : public Sim_Method {
public:
  using rng_t = wcs::RNGen<std::uniform_real_distribution, double>;
  using v_desc_t = Sim_Method::v_desc_t;
  using priority_t = std::pair<reaction_rate_t, v_desc_t>;
  using propensisty_list_t = std::vector<priority_t>;


  SSA_Direct();
  ~SSA_Direct() override;
  /// Initialize propensity list
  void init(std::shared_ptr<wcs::Network>& net_ptr,
            const unsigned max_iter,
            const double max_time,
            const unsigned rng_seed) override;

  /// Main loop of SSA
  std::pair<unsigned, sim_time_t> run() override;

  static bool less(const priority_t& v1, const priority_t& v2);

  rng_t& rgen_e();
  rng_t& rgen_t();

protected:
  void build_propensity_list();
  priority_t& choose_reaction();
  sim_time_t get_reaction_time();
  void update_reactions(priority_t& firing,
                        const Sim_Method::affected_reactions_t& affected);

protected:
  /// Cumulative propensity of reactions events
  propensisty_list_t m_propensity;
  rng_t m_rgen_evt; ///< RNG for events
  rng_t m_rgen_tm; ///< RNG for event times
  /// map from vertex descriptor to propensity
  std::unordered_map<v_desc_t, size_t> m_pindices;
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_SIM_METHODS_SSA_DIRECT_HPP__
