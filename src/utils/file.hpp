/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef  __WCS_UTILS_FILE__
#define  __WCS_UTILS_FILE__
#include <string>

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#if defined(WCS_HAS_STD_FILESYSTEM)
#include <filesystem>
#else
#include <boost/filesystem.hpp>
#endif

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

/// Split a path string into three components: parent dir, stem, and extension
void extract_file_component(const std::string path, std::string& parent_dir,
                            std::string& stem, std::string& extension);

/// Returns a new path string that has a stem appended with the given str
std::string append_to_stem(const std::string path, const std::string str);

/**@}*/
} // end of namespace wcs
#endif //  __WCS_UTILS_FILE__
