/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __STATE_IO_CEREAL_HPP__
#define __STATE_IO_CEREAL_HPP__

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#if defined(WCS_HAS_CEREAL)
#include <cereal/archives/binary.hpp>
#include <iostream>
#include "streamvec.hpp"
#include "streambuff.hpp"
#include "traits.hpp" // is_trivially_copyable

namespace wcs {
template <typename T>
constexpr bool is_custom_bin_cerealizable()
{
  return (!std::is_arithmetic<T>::value &&
           std::is_trivially_copyable<T>::value);
}
} // end of namespace wcs

#define ENABLE_CUSTOM_CEREAL(T) \
namespace cereal { \
  inline std::enable_if_t<wcs::is_custom_bin_cerealizable<T>(), void> \
  CEREAL_SAVE_FUNCTION_NAME(BinaryOutputArchive & ar, T const & t) \
  { \
    ar.saveBinary(std::addressof(t), sizeof(t)); \
  } \
  inline std::enable_if_t<wcs::is_custom_bin_cerealizable<T>(), void> \
  CEREAL_LOAD_FUNCTION_NAME(BinaryInputArchive & ar, T & t) \
  { \
    ar.loadBinary(std::addressof(t), sizeof(t)); \
  } \
}

namespace wcs {

template <typename T>
void save_state(const T& state, std::ostream& os)
{
  // Create an output archive with the given stream
  cereal::BinaryOutputArchive oarchive(os);

  oarchive(state); // Write the data to the archive
  // archive goes out of scope,
  // ensuring all contents are flushed to the stream
}

template <typename T>
void load_state(T& state, std::istream& is)
{
  // Create an input archive using the given stream
  cereal::BinaryInputArchive iarchive(is);

  iarchive(state); // Read the data from the archive
}

} // end of namespace wcs
#endif // WCS_HAS_CEREAL

#endif // __STATE_IO_CEREAL_HPP__
