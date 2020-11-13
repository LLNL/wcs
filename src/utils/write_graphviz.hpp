/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __WCS_UTILS_WRITE_GRAPHVIZ_HPP__
#define __WCS_UTILS_WRITE_GRAPHVIZ_HPP__
#include <iostream>
#include <string>
#include <type_traits>
#include "wcs_types.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

template <typename G>
std::ostream& write_graphviz(std::ostream& os, const G& g,
                             partition_id_t pid = unassigned_partition);

template <typename G>
std::ostream& operator<<(std::ostream& os, const std::pair<const G&, partition_id_t>& gp);


template <typename G>
bool write_graphviz(const std::string& out_filename, const G& g,
                    partition_id_t pid = unassigned_partition);

/**@}*/
} // end of namespace wcs

#include "write_graphviz_impl.hpp"
#endif // __WCS_UTILS_WRITE_GRAPHVIZ_HPP__
