/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include "utils/file.hpp"
#include "utils/exception.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

void extract_file_component(const std::string path,
                            std::string& parent_dir,
                            std::string& stem,
                            std::string& ext)
{
#if defined(WCS_HAS_STD_FILESYSTEM)
  const auto fn = std::filesystem::path(path);
  if (!fn.has_stem()) {
    parent_dir.clear();
    stem.clear();
    ext.clear();
    return;
  }
  stem = std::string(fn.stem());

  if (!fn.has_extension()) {
    ext = "";
  } else {
    ext = std::string(fn.extension());
  }

  if (!fn.has_parent_path()) {
    parent_dir = "";
  } else {
    parent_dir = std::string(fn.parent_path()) + '/';
  }
#else
  const auto fn = boost::filesystem::path(path);
  if (!fn.has_stem()) {
    parent_dir.clear();
    stem.clear();
    ext.clear();
    return;
  }
  stem = fn.stem().string();

  if (!fn.has_extension()) {
    ext = "";
  } else {
    ext = fn.extension().string();
  }

  if (!fn.has_parent_path()) {
    parent_dir = "";
  } else {
    parent_dir = fn.parent_path().string() + '/';
  }
#endif
}

std::string append_to_stem(const std::string path, const std::string str)
{
  std::string parent_dir;
  std::string stem;
  std::string ext;
  extract_file_component(path, parent_dir, stem, ext);
  return (parent_dir + stem + str + ext);
}

/**@}*/
} // end of namespace wcs
