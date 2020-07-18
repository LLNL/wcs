/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef  __WCS_UTILS_INPUT_FILETYPE__
#define  __WCS_UTILS_INPUT_FILETYPE__
#include <string>


namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

class input_filetype {
 public:
  enum input_type { _ioerror_=0, _graphml_, _sbml_, _unknown_ };
  input_filetype(const std::string filename);
  input_type detect() const;

 private:
  std::string m_filename;
};

/**@}*/
} // end of namespace wcs
#endif //  __WCS_UTILS_INPUT_FILETYPE__
