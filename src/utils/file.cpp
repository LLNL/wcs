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

bool check_if_file_exists(const std::string filename)
{
#if defined(WCS_HAS_STD_FILESYSTEM)
  return std::filesystem::exists(std::filesystem::path(filename));
#else
  return boost::filesystem::exists(boost::filesystem::path(filename));
#endif
}

/** Return a library name that ends with '.so' with the rest the same as the
 * model file name */
std::string get_libname_from_model(const std::string& model_filename)
{
  std::string dir, stem, ext;
  extract_file_component(model_filename, dir, stem, ext);

  if (ext == ".so") {
    WCS_THROW("Model filename should not have the extension '.so'");
  }
  return dir + stem + ".so";
}

/** Create a default output filename in case that no name is given.
 *  File name is the same as the model name except that the extension is `.out`
 *  and the file is going be under current working directory.
 */
std::string get_default_outname_from_model(const std::string& model_filename)
{
  std::string dir, stem, ext;
  extract_file_component(model_filename, dir, stem, ext);

  if (ext == ".out") {
    WCS_THROW("Model filename should not have the extension '.out'");
  }
  return stem + ".out";
}

/**@}*/
} // end of namespace wcs
