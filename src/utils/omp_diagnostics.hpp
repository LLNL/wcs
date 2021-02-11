/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef WCS_UTILS_OMP_DIAGNOSTICS_HPP
#define WCS_UTILS_OMP_DIAGNOSTICS_HPP

#include <cstdint>
#include <vector>
#include <string>

namespace wcs {

/**
 *  Obtain the processor affinity of the current thread.
 *  This structure needs to be thread local.
 */
struct my_omp_affinity {
  using cpuid_t = uint8_t;
  int m_tid;
  int m_num_threads;
  int m_my_level;
  std::vector<cpuid_t> m_cpus;
  std::vector<int> m_ancestor_id;

  /**
   * When this is called inside of a parallel region, it will gather
   * information on the processor affinity of the calling thread.
   */ 
  void get();
  /**
   * Print out the processor affinity of the thread gathered by `get()`.
   * This call does not need to be inside of a parallel region as it
   * simply prints out the information gathered.
   */ 
  void print() const;
};

std::string get_omp_version();
std::string to_string_omp_schedule_kind(int kind);
void set_static_schedule();

} // namespace wcs
#endif // WCS_UTILS_OMP_DIAGNOSTICS_HPP
