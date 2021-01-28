/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef WCS_WCS_PARAMS_HPP
#define WCS_WCS_PARAMS_HPP

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
//#error "no config"
#endif

namespace wcs {

struct cmd_line_opts {
  std::string m_input_model;
  std::string m_all_setup;
  std::string m_sim_setup;
  std::string m_part_setup;
  std::string m_des_setup;
  bool m_is_set;

  bool parse_cmd_line(int argc, char** argv);
  void show() const;
};

} // end of namespace wcs

#endif // WCS_WCS_PARAMS_HPP
