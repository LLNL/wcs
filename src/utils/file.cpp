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
#include <sys/stat.h> // mkdir(), open()
#include <cerrno> // errno
#include <unistd.h> // getuid(), fsync(), close()
#include <cstdio> // perror()
#include <libgen.h> // dirname()
#include <cstring> // strncpy()
#include <fcntl.h> // O_RDONLY
#include <fstream> // filebuf

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

/**
 * For a base path "/a/b/c" and a full path "/a/b/c/d/f.cpp", return "d/f.cpp".
 * This is useful if the current working directory is the base path.
 */
std::string get_subpath(const std::string& basepath,
                        const std::string& fullpath)
{
  std::string dir, stem, ext;
  extract_file_component(fullpath, dir, stem, ext);
  auto pos = dir.find(basepath);

  if  (pos == std::string::npos) {
    return fullpath;
  }

  auto dir_match = dir.substr(pos);

  if (dir_match.length() < basepath.length()) {
    return fullpath;
  }
  if (((dir_match.length() == basepath.length()) && (basepath.back() == '/')) ||
      ((dir_match.length() == basepath.length()+1) && (dir_match.back() == '/')))
  {
    return stem + ext;
  }

  return dir_match.substr(basepath.length()) + stem + ext;
}

/**
 * Create a directory `path` with the given access permission `m`.
 * It is the same as the POSIX mkdir() except that it does not
 * return as failure if the directory to create already exists with
 * the same access permission. However, the access group name is not
 * checked. The default mode of creation is 0700.
 */
int mkdir_as_needed (const std::string& path, const mode_t m)
{
  if (path.empty()) {
    errno = ENOENT;
    perror(path.c_str());
    return -1;
  }

  struct stat sb;
  const mode_t RWX_UGO = (S_IRWXU | S_IRWXG | S_IRWXO);

  if ((stat(path.c_str(), &sb) == 0) &&  // already exists
      (S_ISDIR(sb.st_mode) != 0) &&      // is a directory
      (sb.st_uid == getuid()) &&         // same owner
      (((sb.st_mode & RWX_UGO) ^ (m & RWX_UGO)) ==
       static_cast<mode_t>(00))) // same permission
  {
    // Assuming the access group are the same
    return 0; // report success as the same already exists
  }
  std::cerr << "Creating a directory '" + path + "'" << std::endl;

  const mode_t old_mask = umask (0);
  if (mkdir (path.c_str(), m) != 0) {
    perror(path.c_str());
    return -1;
  }
  umask (old_mask);

  return 0;
}

bool sync_directory(const std::string& path)
{ // Flush new directory entry https://lwn.net/Articles/457671/
  char path_copy [PATH_MAX+1] = {'\0'};
  int odir_fd = -1;
  char *odir = NULL;
  bool rc = true;

  strncpy(path_copy, path.c_str(), PATH_MAX);
  odir = dirname(path_copy);

  if ((odir_fd = open(odir, O_RDONLY)) < 0) {
    std::cerr << "Cannot open the directory " + std::string(odir) << std::endl;
    rc = false;
  } else {
    if (fsync(odir_fd) < 0) {
      std::cerr << "Cannot flush the directory " + std::string(odir) << std::endl;
      rc = false;
    }
    if (close(odir_fd) < 0) {
      std::cerr << "Cannot close the directory " + std::string(odir) << std::endl;
      rc = false;
    }
  }
  return rc;
}

// https://stackoverflow.com/questions/676787/how-to-do-fsync-on-an-ofstream
void fsync_ofstream(std::ofstream& os)
{
  class my_filebuf : public std::filebuf
  {
   public:
    int handle() { return _M_file.fd(); }
  };

  os.flush();
  fsync(static_cast<my_filebuf&>(*os.rdbuf()).handle());
}

/**@}*/
} // end of namespace wcs
