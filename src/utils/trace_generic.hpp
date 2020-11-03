/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef	 __WCS_UTILS_TRACE_GENERIC_HPP__
#define	 __WCS_UTILS_TRACE_GENERIC_HPP__
#include <list>
#include "wcs_types.hpp"
#include "utils/trajectory.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

class TraceGeneric : public Trajectory
{
public:
  using tentry_t = typename std::pair<sim_time_t, cnt_updates_t>;
  using trace_t = typename std::list<tentry_t>;

  TraceGeneric(const std::shared_ptr<wcs::Network>& net_ptr);
  TraceGeneric(const TraceGeneric& other) = default;
  TraceGeneric(TraceGeneric&& other) = default;
  TraceGeneric& operator=(const TraceGeneric& other) = default;
  TraceGeneric& operator=(TraceGeneric&& other) = default;

  ~TraceGeneric() override;
  void record_step(const sim_time_t t, cnt_updates_t&& updates) override;
  void finalize() override;

protected:
  size_t estimate_tmpstr_size() const;
  std::ostream& write_header(std::ostream& os) const override;
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
#endif // __WCS_UTILS_TRACE_GENERIC_HPP__
