/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef _WCS_UTILS_TIMER_HPP_
#define _WCS_UTILS_TIMER_HPP_

#include <chrono>

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

inline double get_time() {
  using namespace std::chrono;
  return duration_cast<duration<double>>(
           steady_clock::now().time_since_epoch()).count();
}

/**@}*/
} // namespace wcs

#endif  // _WCS_UTILS_TIMER_HPP_
