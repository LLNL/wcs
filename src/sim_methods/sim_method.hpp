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
#include <cmath>
#include <limits>
#include <unordered_map>
#include "wcs_types.hpp"
#include "reaction_network/network.hpp"
#include "utils/rngen.hpp"
#include "utils/trace_ssa.hpp"
#include "utils/samples.hpp"

namespace wcs {
/** \addtogroup wcs_sim_methods
 *  *  @{ */

class Sim_Method {
public:
  using v_desc_t = wcs::Network::v_desc_t;
  using trace_t = wcs::TraceSSA;
  using samples_t = wcs::Samples;
  using sim_time_t = wcs::sim_time_t;
  using reaction_rate_t = wcs::reaction_rate_t;

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
  /// Record the final state of simulation
  void record_final_state(const sim_time_t dt, const v_desc_t rv);

  /// Check whether to record the state at current step
  bool check_to_record();
  /// Record the state at current step as needed
  bool check_to_record(const sim_time_t dt, const v_desc_t rv);

  virtual std::pair<sim_iter_t, sim_time_t> run() = 0;

  /// Allow access to the internal tracer
  trace_t& trace();
  /// Allow access to the internal sampler
  samples_t& samples();

protected:

  /** The pointer to the reaction network being monitored.
   *  Make sure the network object does not get destroyed
   *  while the trace refers to it.
   */
  std::shared_ptr<wcs::Network> m_net_ptr;

  sim_iter_t m_max_iter; ///< Upper bound on simulation iteration
  sim_time_t m_max_time; ///< Upper bound on simulation time

  sim_iter_t m_cur_iter; ///< Current simulation iteration
  sim_time_t m_sim_time; ///< Current simulation time

  bool m_enable_tracing; ///< Whether to enable tracing
  bool m_enable_sampling; ///< Whether to enable sampling

  sim_iter_t m_sample_iter_interval;
  sim_time_t m_sample_time_interval;

  sim_iter_t m_next_sample_iter; ///< Next iteration to sample
  sim_time_t m_next_sample_time; ///< Next time to sample

  trace_t m_trace; ///< Tracing record
  samples_t m_samples; ///< Sampling record

  sim_time_t dt_sample; ///< time elapsed since the last sample
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_SIM_METHODS_SIM_METHOD_HPP__
