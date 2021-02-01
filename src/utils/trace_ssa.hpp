/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef	 __WCS_UTILS_TRACE_SSA_HPP__
#define	 __WCS_UTILS_TRACE_SSA_HPP__
#include <list>
#include <vector>
#include "utils/trajectory.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

class TraceSSA : public Trajectory {
public:
  using tentry_t = typename std::pair<sim_time_t, r_desc_t>;
  using trace_t = typename std::list<tentry_t>;

  TraceSSA(const std::shared_ptr<wcs::Network>& net_ptr);
  TraceSSA(const TraceSSA& other) = default;
  TraceSSA(TraceSSA&& other) = default;
  TraceSSA& operator=(const TraceSSA& other) = default;
  TraceSSA& operator=(TraceSSA&& other) = default;

  ~TraceSSA() override;
  void initialize() override;
  using Trajectory::record_step;
  void record_step(const sim_time_t t, const r_desc_t r) override;
  void finalize(const sim_time_t t) override;

protected:
  std::ostream& write_header(std::ostream& os) const override;
  size_t estimate_tmpstr_size() const;
  std::ostream& print_stats(const sim_time_t sim_time,
                            const std::string rlabel,
                            std::string& tmpstr, std::ostream& os) const;
  std::ostream& write(std::ostream& os) override;
  void flush() override;

protected:
  /// Trace records
  trace_t m_trace;
  /// Show how many times each reaction fires
  std::vector<r_cnt_t> m_reaction_counts;
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_UTILS_TRACE_SSA_HPP__
