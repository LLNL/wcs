#ifndef _WCS_UTILS_TIMER_HPP_
#define _WCS_UTILS_TIMER_HPP_

#include <chrono>

namespace wcs {

inline double get_time() {
  using namespace std::chrono;
  return duration_cast<duration<double>>(
           steady_clock::now().time_since_epoch()).count();
}

} // namespace wcs

#endif  // _WCS_UTILS_TIMER_HPP_
