/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef WCS_SSA_PARAMS_HPP
#define WCS_SSA_PARAMS_HPP

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#include <string>
#include "wcs_types.hpp"

namespace wcs {
/** \addtogroup wcs_params
 *  @{ */

struct SSA_Params {
  SSA_Params();
  void getopt(int& argc, char** &argv);
  void print_usage(const std::string exec, int code);
  void print() const;

  unsigned m_seed;
  wcs::sim_iter_t m_max_iter;
  wcs::sim_time_t m_max_time;
  int m_method;
  bool m_tracing;
  bool m_sampling;
  wcs::sim_iter_t m_iter_interval;
  wcs::sim_time_t m_time_interval;
  unsigned m_frag_size;
  bool m_is_frag_size_set;

  std::string m_infile;
  std::string m_outfile;
  std::string m_gvizfile;

  bool m_is_iter_set;
  bool m_is_time_set;
};

/**@}*/
} // end of namespace wcs
#endif // WCS_SSA_PARAMS_HPP
