/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __WCS_SIM_METHODS_SIM_METHOD_HPP__
#define __WCS_SIM_METHODS_SIM_METHOD_HPP__

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#include <cmath>
#include <limits>
#include <unordered_map>
#include "wcs_types.hpp"
#include "reaction_network/network.hpp"
#include "utils/rngen.hpp"
#include "utils/trace_ssa.hpp"
#include "utils/samples.hpp"
#include "sim_methods/sim_state_change.hpp"

namespace wcs {
/** \addtogroup wcs_sim_methods
 *  @{ */

class Sim_Method {
public:
  using v_desc_t = wcs::Network::v_desc_t;
  using trace_t = wcs::TraceSSA;
  using samples_t = wcs::Samples;
  using sim_time_t = wcs::sim_time_t;
  using reaction_rate_t = wcs::reaction_rate_t;

  /** Type for keeping track of species updates to facilitate undoing
   *  reaction processing.  */
  using update_t = Sim_State_Change::update_t;
  using update_list_t = Sim_State_Change::update_list_t;
  /** Type for the list of reactions that share any of the species with the
   *  firing reaction */
  using affected_reactions_t = Sim_State_Change::affected_reactions_t;

  enum result_t {Success, Empty, Inactive};

  Sim_Method();
  virtual ~Sim_Method();
  virtual void init(std::shared_ptr<wcs::Network>& net_ptr,
                    const sim_iter_t max_iter,
                    const double max_time,
                    const unsigned rng_seed) = 0;

  /// Enable tracing to record state at every event
  void set_tracing();
  /// Disable tracing to record state at every event
  void unset_tracing();
  /// Enable sampling to record state at every given time interval
  void set_sampling(const sim_time_t time_interval);
  /// Enable sampling to record state at every given iteration interval
  void set_sampling(const sim_iter_t iter_interval);
  /// Disable sampling
  void unset_sampling();

  /// Record the initial state of simulation
  void record_initial_state(const std::shared_ptr<wcs::Network>& net_ptr);

  /// Record the state at current step
  void record(const v_desc_t rv);
  /// Record the state at time t. This allows tracing/sampling using history.
  void record(const sim_time_t t, const v_desc_t rv);

 #if defined(WCS_HAS_ROSS)
  /**
   * Record as many states as the given number of iterations from the beginning
   * of the digest list */
  virtual void record_first_n(const sim_iter_t num) = 0;
 #endif // defined(WCS_HAS_ROSS)

  virtual std::pair<sim_iter_t, sim_time_t> run() = 0;

  /// Allow access to the internal tracer
  trace_t& trace();
  /// Allow access to the internal sampler
  samples_t& samples();

  bool fire_reaction(Sim_State_Change& digest);

  void undo_species_updates(const update_list_t& updates) const;
  bool undo_reaction(const Sim_Method::v_desc_t& rd_undo) const;

protected:

  /** The pointer to the reaction network being monitored.
   *  Make sure the network object does not get destroyed
   *  while the trace refers to it.
   */
  std::shared_ptr<wcs::Network> m_net_ptr;

  sim_iter_t m_max_iter; ///< Upper bound on simulation iteration
  sim_time_t m_max_time; ///< Upper bound on simulation time

  sim_iter_t m_sim_iter; ///< Current simulation iteration
  sim_time_t m_sim_time; ///< Current simulation time

  bool m_enable_tracing; ///< Whether to enable tracing
  bool m_enable_sampling; ///< Whether to enable sampling

  trace_t m_trace; ///< Tracing record
  samples_t m_samples; ///< Sampling record
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_SIM_METHODS_SIM_METHOD_HPP__
