/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef PRINT_VERTICES_HPP
#define PRINT_VERTICES_HPP
#include <string>
#include "reaction_network/network.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

void print_vertices(const std::shared_ptr<wcs::Network>& net_ptr,
                    const std::string& outfile_name);

/**@}*/
} // end of namespasce wcs
#endif // PRINT_VERTICES_HPP

