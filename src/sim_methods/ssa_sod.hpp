/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __WCS_SIM_METHODS_SSA_SOD_HPP__
#define __WCS_SIM_METHODS_SSA_SOD_HPP__

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#include <cmath>
#include <limits>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/identity.hpp>
#include "sim_methods/sim_method.hpp"

namespace wcs {
/** \addtogroup wcs_sim_methods
 *  @{ */

/**
 * Sorted optimized direct method.
 * This method takes advantage of a dependency graph similarly to the optimized
 * direct method, and it maintains the propensity list in sorted order.
 */
class SSA_SOD : public Sim_Method {
public:
  using rng_t = wcs::RNGen<std::uniform_real_distribution, double>;
  using v_desc_t = Sim_Method::v_desc_t;

  struct priority_t {
    reaction_rate_t m_rate; ///< reaction rate
    /**
     * Cumulative reaction rate.
     * It is defined as a mutable field such that it can be directly modified
     * without the overhead of using the modifier which would evaluate whether
     * to reposition the item after a change. Skipping it is only acceptable
     * as the key comparison for indexing does not rely on the field and, thus,
     * no repositioning occurs by updating its value.
     */
    mutable reaction_rate_t m_curate;
    v_desc_t m_rvd; ///< reaction vertex descriptor

    /**
     *  TODO: If no operator'<' is defined for the vertex descriptor type,
     *  it should rely on the comparison of the hash values of the operands
     *  when possible.
     */
    bool operator<(const priority_t& rhs) const {
      return (m_rate < rhs.m_rate) ||
             ((m_rate == rhs.m_rate) && (m_rvd > rhs.m_rvd));
    }
  };

  struct tag_rate {};
  struct tag_rvd {};

  using propensity_list_t = boost::multi_index_container<
    priority_t,
    boost::multi_index::indexed_by<
      boost::multi_index::ordered_non_unique<
        boost::multi_index::tag<tag_rate>,
        boost::multi_index::identity<priority_t>
      >,
      boost::multi_index::hashed_unique<
        boost::multi_index::tag<tag_rvd>,
        boost::multi_index::member<priority_t, v_desc_t, &priority_t::m_rvd>
      >
    >
  >;

  using rate_idx_t = propensity_list_t::index<tag_rate>::type;
  using rvd_idx_t = propensity_list_t::index<tag_rvd>::type;


  SSA_SOD(const std::shared_ptr<wcs::Network>& net_ptr);
  SSA_SOD(SSA_SOD&& other) = default;
  SSA_SOD& operator=(SSA_SOD&& other) = default;
  ~SSA_SOD() override;

  /// Initialize propensity list
  void init(const unsigned max_iter,
            const double max_time,
            const unsigned rng_seed) override;

  /**
   * Determines when the next reaction to occur.
   * When successful, this function returns Sim_Method::Success. Otherwise,
   * it returns a failure code.
   */
  Sim_Method::result_t schedule(sim_time_t& t);
  /**
   * Determine which reaction to fire and execute it at the given time.
   * Check the simulation termination condition at the beginning. If it is not
   * to be terminated yet, proceed and return true. Otherwise, stop immediately
   * and return false.
   */
  bool forward(const sim_time_t t);
  /// Main loop of SSA
  std::pair<unsigned, sim_time_t> run() override;

 #if defined(WCS_HAS_ROSS)
  void backward(Sim_State_Change& digest);

  /** Record as many states as the given number of iterations from the
   *  beginning of the digest list */
  void record_first_n(const sim_iter_t num) override;
 #endif // defined(WCS_HAS_ROSS)

  rng_t& rgen_e();
  rng_t& rgen_t();

protected:
  void build_propensity_list();
  priority_t choose_reaction();
  sim_time_t get_reaction_time();
  void update_reactions(const v_desc_t& rd_fired,
                        const Sim_Method::affected_reactions_t& affected,
                        const bool check_reaction);

  void save_rgen_state(Sim_State_Change& digest) const;
  void load_rgen_state(const Sim_State_Change& digest);

protected:
  /// Cumulative propensity of reactions events
  propensity_list_t m_propensity;
  rng_t m_rgen_evt; ///< RNG for events
  rng_t m_rgen_tm; ///< RNG for event times

 #if defined(WCS_HAS_ROSS)
  using digest_list_t = std::list<Sim_State_Change>;
  digest_list_t m_digests;
 #endif // defined(WCS_HAS_ROSS)
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_SIM_METHODS_SSA_SOD_HPP__
