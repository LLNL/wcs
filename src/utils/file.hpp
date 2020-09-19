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
#include <filesystem>
// TODO: switch to boost:filesystem if filesystem is not supported
// C++ preprocessor definition set by cmake if cmake tests true

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

void extract_file_component(const std::string path, std::string& parent_dir,
                            std::string& stem, std::string& extension);

/**@}*/
} // end of namespace wcs
#endif //  __WCS_UTILS_FILE__
