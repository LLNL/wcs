/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __WCS_REACTION_NETWORK_RATE_STATS_HPP__
#define __WCS_REACTION_NETWORK_RATE_STATS_HPP__

#include <limits>
#include <string>
#include <iostream>
#include "wcs_types.hpp"
#include "utils/to_string.hpp"

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  @{ */

struct RateStats {
  using rrate_t = reaction_rate_t;
  rrate_t m_rmin;
  rrate_t m_rmax;
  rrate_t m_rsum;
  double m_ravg;
  size_t m_imin;
  size_t m_imax;
  size_t m_isum;
  double m_iavg;
  size_t m_num_reactions;

  void init();
  void read(const rrate_t val, const size_t n_inputs);
  void set_avg();
  void print() const;
};

inline void RateStats::init()
{
  m_rmin = std::numeric_limits<rrate_t>::max();
  m_rmax = std::numeric_limits<rrate_t>::min();
  m_rsum = static_cast<rrate_t>(0);
  m_ravg = 0.0;
  m_imin = std::numeric_limits<size_t>::max();
  m_imax = std::numeric_limits<size_t>::min();
  m_isum = 0ul;
  m_iavg = 0.0;
  m_num_reactions = 0ul;
}

inline void RateStats::read(const rrate_t val, const size_t num_inputs)
{
  if (val < m_rmin) m_rmin = val;
  if (val > m_rmax) m_rmax = val;
  m_rsum += val;

  if (num_inputs < m_imin) m_imin = num_inputs;
  if (num_inputs > m_imax) m_imax = num_inputs;
  m_isum += num_inputs;

  m_num_reactions ++;
}

inline void RateStats::set_avg()
{
  m_ravg = m_rsum / static_cast<double>(m_num_reactions);
  m_iavg = m_isum / static_cast<double>(m_num_reactions);
}

inline void RateStats::print() const
{
  if (m_num_reactions == 0ul) {
    return;
  }
  std::string msg = "Reaction statistic:";
  msg += "\n - min of rate: " + to_string_in_scientific(m_rmin);
  msg += "\n - max of rate: " + to_string_in_scientific(m_rmax);
  msg += "\n - sum of rate: " + to_string_in_scientific(m_rsum);
  msg += "\n - avg of rate: " + to_string_in_scientific(m_ravg);
  msg += "\n - min of number of inputs: " + std::to_string(m_imin);
  msg += "\n - max of number of inputs: " + std::to_string(m_imax);
  msg += "\n - sum of number of inputs: " + std::to_string(m_isum);
  msg += "\n - avg of number of inputs: " + to_string_in_scientific(m_iavg);
  msg += "\n - number of reactions: " + std::to_string(m_num_reactions);
  msg += "\n";

 #if defined(_OPENMP)
  #pragma omp critical
 #endif // defined(_OPENMP)
  {
    std::cout << msg;
  }
}

/**@}*/
} // end of namespace wcs
#endif // __WCS_REACTION_NETWORK_RATE_STATSE_HPP__
