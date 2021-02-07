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

  unsigned seed;
  wcs::sim_iter_t max_iter;
  wcs::sim_time_t max_time;
  int method;
  bool tracing;
  bool sampling;
  wcs::sim_iter_t iter_interval;
  wcs::sim_time_t time_interval;
  unsigned frag_size;
  bool is_frag_size_set;

  std::string infile;
  std::string outfile;
  std::string gvizfile;

  bool is_iter_set;
  bool is_time_set;
};

/**@}*/
} // end of namespace wcs
#endif // WCS_SSA_PARAMS_HPP
