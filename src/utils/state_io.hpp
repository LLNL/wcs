/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __STATE_IO_HPP__
#define __STATE_IO_HPP__
#include <iostream>
#include "streamvec.hpp"
#include "streambuff.hpp"
#include "traits.hpp"

namespace wcs {

/**
 * Override the stream operators only for trivially copyable types.
 * Users call `bits()` on the object of such a type to use the overriden
 * interfaces. Note that the function is only defined for trivially
 * copyable ones. Calling on an object of a wrong type would generate a
 * compiler error.
 * https://stackoverflow.com/questions/1559254/are-there-binary-memory-streams-in-c
 */
template<typename T> struct bits_t {
  using value_type = typename std::remove_reference<T>::type;
  T v;
};

/// Return wrapper struct bits_t of a trivially copyable type
template <typename T>
typename std::enable_if<
  !is_vector<T>::value &&
  std::is_trivially_copyable<T>::value,
  bits_t<T&> >::type
bits(T &v)
{
  return bits_t<T&>{v};
}

/// Return wrapper struct bits_t of a trivially copyable read-only type
template <typename T>
typename std::enable_if<
  !is_vector<T>::value &&
  std::is_trivially_copyable<T>::value,
  bits_t<const T&> >::type
bits(const T& v)
{
  return bits_t<const T&>{v};
}

/// Return wrapper struct bits_t of the vector of a trivially copyable type
template <typename T>
typename std::enable_if<
  is_vector<T>::value &&
  !is_bool<typename T::value_type>::value &&
  std::is_trivially_copyable<typename T::value_type>::value,
  bits_t<T&> >::type
bits(T &v)
{
  return bits_t<T&>{v};
}

/** Return wrapper struct bits_t of the vector of trivially copyable read-only
 *  type */
template <typename T>
inline typename std::enable_if<
  is_vector<T>::value &&
  !is_bool<typename T::value_type>::value &&
  std::is_trivially_copyable<typename T::value_type>::value,
  bits_t<const T&> >::type
bits(const T& v)
{
  return bits_t<const T&>{v};
}

/// Output data wrapped in bits_t struct into the given stream
template<typename S, typename T>
S& operator<<(S &os, const bits_t<T&>& b);

/// Load the data from the given stream into data space wrapped in bits_t struct
template<typename S, typename T>
S& operator>>(S& is, const bits_t<T&>& b);


template<typename ObjT,
         typename CharT = char,
         typename Traits = std::char_traits<CharT> >
bool save_state(const ObjT& obj, std::vector<CharT>& buffer);

template<typename ObjT,
         typename CharT = char,
         typename Traits = std::char_traits<CharT> >
bool load_state(ObjT& obj, const std::vector<CharT>& buffer);

template<typename ObjT,
         typename CharT = char,
         typename Traits = std::char_traits<CharT> >
bool save_state(const ObjT& obj, CharT* buffer);

template<typename ObjT,
         typename CharT = char,
         typename Traits = std::char_traits<CharT> >
bool load_state(ObjT& obj, const CharT* buffer);

} // end of namespace wcs

#include "state_io_impl.hpp"

#endif // __STATE_IO_HPP__
