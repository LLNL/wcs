/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef	 __WCS_UTILS_TRACE_ODE_HPP__
#define	 __WCS_UTILS_TRACE_ODE_HPP__
#include <vector>
#include "wcs_types.hpp"
#include "utils/trajectory.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

class TraceODE : Trajectory
{
public:
  using tentry_t = typename std::pair<sim_time_t, update_list_t>;
  using trace_t = typename std::list<tentry_t>;

  ~TraceODE() override;
  void record_step(const sim_time_t t, const update_list_t& updates) override;
  void finalize() override;

protected:
  void build_index_maps() override;
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
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_UTILS_TRACE_ODE_HPP__
