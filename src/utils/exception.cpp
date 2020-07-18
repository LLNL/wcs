/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include "utils/exception.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

exception::exception(const std::string msg)
: m_message(msg)
{}

const char* exception::what() const noexcept
{
  return m_message.c_str();
}

std::ostream& operator<<(std::ostream& os, const exception& e)
{
  using std::operator<<;
  os << e.what();
  return os;
}

/**@}*/
} // end of namespace wcs
