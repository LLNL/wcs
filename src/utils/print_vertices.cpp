/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include <iostream>
#include <fstream>
#include "utils/exception.hpp"
#include "utils/print_vertices.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

void print_vertices(const std::shared_ptr<wcs::Network>& net_ptr,
                    const std::string& outfile_name)
{
  if (outfile_name.empty()) {
    return;
  }

  if (!net_ptr) {
    WCS_THROW("Invaid pointer for reaction network.");
  }
  const wcs::Network::graph_t& g = net_ptr->graph();

  {
    size_t idx = 0ul;
    
    std::ofstream os(outfile_name + ".species_index.txt");
    using std::operator<<;
    
    std::string str;
    str.reserve(2048);
    for (const auto& sd : net_ptr->species_list()) {
      str = std::to_string(idx ++) + ' ' + g[sd].get_label();
      os << str << std::endl;
    }
  }

  {
    size_t idx = 0ul;
    
    std::ofstream os(outfile_name + ".reaction_index.txt");
    using std::operator<<;
    
    std::string str;
    str.reserve(2048);
    for (const auto& rd : net_ptr->reaction_list()) {
      str = std::to_string(idx ++) + ' ' + g[rd].get_label();
      os << str << std::endl;
    }
  }
}

/**@}*/
} // end of namespace wcs
