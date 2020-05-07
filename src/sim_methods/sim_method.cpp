/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include "sim_methods/sim_method.hpp"
#include <limits>

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  *  @{ */

Sim_Method::Sim_Method()
: m_net_ptr(nullptr),
  m_max_iter(static_cast<sim_iter_t>(0u)),
  m_max_time(static_cast<sim_time_t>(0)),
  m_cur_iter(static_cast<sim_iter_t>(0u)),
  m_sim_time(static_cast<sim_time_t>(0)),
  m_enable_tracing(false),
  m_enable_sampling(false),
  m_sample_iter_interval(static_cast<sim_iter_t>(0u)),
  m_sample_time_interval(static_cast<sim_time_t>(0)),
  m_next_sample_iter(static_cast<sim_iter_t>(0u)),
  m_next_sample_time(static_cast<sim_time_t>(0)),
  dt_sample(static_cast<sim_time_t>(0))
{
  using directed_category
    = typename boost::graph_traits<wcs::Network::graph_t>::directed_category;

  constexpr bool is_bidirectional
    = std::is_same<directed_category, boost::bidirectional_tag>::value;

  if constexpr (!is_bidirectional) {
    WCS_THROW("Cannot get species population without in-edges.");
  }
}

Sim_Method::~Sim_Method() {}

void Sim_Method::set_tracing()
{
  m_enable_tracing = true;
  m_enable_sampling = false;
}

void Sim_Method::unset_tracing()
{
  m_enable_tracing = false;
}

void Sim_Method::set_sampling(const sim_time_t time_interval)
{
  m_enable_tracing = false;
  m_enable_sampling = true;
  m_sample_time_interval = time_interval;
  m_next_sample_time = m_sim_time + time_interval;
  m_next_sample_iter = std::numeric_limits<sim_iter_t>::max();
}

void Sim_Method::set_sampling(const sim_iter_t iter_interval)
{
  m_enable_tracing = false;
  m_enable_sampling = true;
  m_sample_iter_interval = iter_interval;
  m_next_sample_iter = m_cur_iter + iter_interval;
  m_next_sample_time = std::numeric_limits<sim_time_t>::infinity();
}

void Sim_Method::unset_sampling()
{
  m_enable_sampling = false;
}

void Sim_Method::record_initial_state(const std::shared_ptr<wcs::Network>& net_ptr)
{ // record initial state of the network
  if (m_enable_tracing) {
    m_trace.record_initial_condition(m_net_ptr);
  } else if (m_enable_sampling) {
    m_samples.record_initial_condition(m_net_ptr);
  }
}

void Sim_Method::record_final_state(const sim_time_t dt, const v_desc_t rv)
{
  if (m_enable_tracing) {
    m_trace.record_reaction(dt, rv);
  } else if (m_enable_sampling) {
    m_samples.record_reaction(rv);
    dt_sample += dt;
    m_samples.take_sample(dt_sample);
  }
}

bool Sim_Method::check_to_record(const sim_time_t dt, const v_desc_t rv)
{
  if (m_enable_tracing) {
    m_trace.record_reaction(dt, rv);
    return true;
  } else if (m_enable_sampling) {
    m_samples.record_reaction(rv);
    dt_sample += dt;
    if (m_cur_iter >= m_next_sample_iter) {
      m_next_sample_iter += m_sample_iter_interval;
      m_samples.take_sample(dt_sample);
      dt_sample = 0.0;
      return true;
    } else if (m_sim_time >= m_next_sample_time) {
      m_next_sample_time += m_sample_time_interval;
      m_samples.take_sample(dt_sample);
      dt_sample = 0.0;
      return true;
    } else {
      return false;
    }
  }
  return false;
}

bool Sim_Method::check_to_record()
{
  if (m_enable_tracing) {
    return true;
  } else if (m_enable_sampling) {
    if (m_cur_iter >= m_next_sample_iter) {
      m_next_sample_iter += m_sample_iter_interval;
      return true;
    } else if (m_sim_time >= m_next_sample_time) {
      m_next_sample_time += m_sample_time_interval;
      return true;
    } else {
      return false;
    }
  }
  return false;
}

Sim_Method::trace_t& Sim_Method::trace() {
  return m_trace;
}

Sim_Method::samples_t& Sim_Method::samples() {
  return m_samples;
}


/**@}*/
} // end of namespace wcs
