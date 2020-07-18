/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef  __WCS_UTILS_EXCEPTION__
#define  __WCS_UTILS_EXCEPTION__
#include <string>
#include <iostream>
#include <exception>

#define WCS_THROW(_MSG_)                                   \
  throw wcs::exception(std::string( __FILE__) + " : line " \
                       + std::to_string(__LINE__) + " : "  \
                       + _MSG_ + '\n')

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

class exception : public std::exception {
 public:
  exception(const std::string message = "");
  const char* what() const noexcept override;

 private:
  std::string m_message;
};

using exception = ::wcs::exception;

std::ostream& operator<<(std::ostream& os, const exception& e);

/**@}*/
} // end of namespace wcs
#endif //  __WCS_UTILS_EXCEPTION__
