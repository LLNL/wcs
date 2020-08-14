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


  SSA_SOD();
  ~SSA_SOD() override;
  /// Initialize propensity list
  void init(std::shared_ptr<wcs::Network>& net_ptr,
            const unsigned max_iter,
            const double max_time,
            const unsigned rng_seed) override;

  /// Main loop of SSA
  std::pair<unsigned, sim_time_t> run() override;

  rng_t& rgen_e();
  rng_t& rgen_t();

protected:
  void build_propensity_list();
  priority_t choose_reaction();
  sim_time_t get_reaction_time();
  void update_reactions(const priority_t& fired,
                        const Sim_Method::affected_reactions_t& affected);

protected:
  /// Cumulative propensity of reactions events
  propensity_list_t m_propensity;
  rng_t m_rgen_evt; ///< RNG for events
  rng_t m_rgen_tm; ///< RNG for event times
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_SIM_METHODS_SSA_SOD_HPP__
