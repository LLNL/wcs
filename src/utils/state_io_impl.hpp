/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __STATE_IO_IMPL_HPP__
#define __STATE_IO_IMPL_HPP__

namespace wcs {

template<typename S, typename T>
S& operator<<(S &os, const bits_t<T&>& b)
{
  if constexpr (is_vector<T>::value) {
    const typename T::size_type sz = b.v.size();
    os.write(reinterpret_cast<const typename S::char_type*>(&sz), sizeof(sz));
    if (sz > 0ul) {
      os.write(reinterpret_cast<const typename S::char_type*>(b.v.data()),
               sz*sizeof(typename T::value_type));
    }
  } else {
    os.write(reinterpret_cast<const typename S::char_type*>(&b.v), sizeof(T));
  }
  return os;
}

template<typename S, typename T>
S& operator>>(S& is, const bits_t<T&>& b)
{
  if constexpr (is_vector<T>::value) {
    typename T::size_type sz = 0ul;
    is.read(reinterpret_cast<typename S::char_type*>(&sz), sizeof(sz));
    b.v.resize(sz);
    if (sz > 0ul) {
      is.read(reinterpret_cast<typename S::char_type*>(b.v.data()),
              sz*sizeof(typename T::value_type));
    }
  } else {
    is.read(reinterpret_cast<typename S::char_type*>(&b.v), sizeof(T));
  }
  return is;
}

template<typename ObjT,
         typename CharT,
         typename Traits>
bool save_state(const ObjT& obj, std::vector<CharT>& buffer) {
  static_assert(std::is_arithmetic<CharT>::value &&
                !std::is_same<CharT, bool>::value,
                "Invalid vector element type");
  /* Resize or reserve vector space to avoid overhead of reallocation
     Especially when there are multiple items to pack and the sizes are
     known in advance. For ostreamvec, the vector is considered empty
     even if the size is larger than 0.
  if (buffer.size()*sizeof(CharT) < sizeof(obj)) {
    buffer.resize((sizeof(obj) + sizeof(CharT) - 1u)/sizeof(CharT));
  }
  */

  ostreamvec<CharT, Traits> ostrmbuf(buffer);
  //streamvec<CharT, Traits> ostrmbuf(buffer);
  std::basic_ostream<CharT, Traits> oss(&ostrmbuf);

  oss << bits(obj);

  return oss.good();
}

template<typename ObjT,
         typename CharT,
         typename Traits>
bool load_state(ObjT& obj, const std::vector<CharT>& buffer) {
  static_assert(std::is_arithmetic<CharT>::value &&
                !std::is_same<CharT, bool>::value,
                "Invalid vector element type");
  istreamvec<CharT, Traits> istrmbuf(buffer);
  //streamvec<CharT, Traits> istrmbuf(buffer, true);
  std::basic_istream<CharT, Traits> iss(&istrmbuf);

  if (istrmbuf.size()*sizeof(CharT) < sizeof(ObjT)) {
    return false;
  }

  iss >> bits(obj);

  return iss.good();
}

template<typename ObjT,
         typename CharT,
         typename Traits>
bool save_state(const ObjT& obj, CharT* buffer) {
  static_assert(std::is_arithmetic<CharT>::value &&
                !std::is_same<CharT, bool>::value,
                "Invalid vector element type");
  ostreambuff<CharT, Traits> ostrmbuf(buffer, sizeof(ObjT));
  //streambuff<CharT, Traits> ostrmbuf(buffer);
  std::basic_ostream<CharT, Traits> oss(&ostrmbuf);

  oss << bits(obj);

  return oss.good();
}

template<typename ObjT,
         typename CharT,
         typename Traits>
bool load_state(ObjT& obj, const CharT* buffer) {
  static_assert(std::is_arithmetic<CharT>::value &&
                !std::is_same<CharT, bool>::value,
                "Invalid vector element type");
  istreambuff<CharT, Traits> istrmbuf(buffer, sizeof(ObjT));
  //streambuff<CharT, Traits> istrmbuf(buffer, true);
  std::basic_istream<CharT, Traits> iss(&istrmbuf);

  if (istrmbuf.size()*sizeof(CharT) < sizeof(ObjT)) {
    return false;
  }

  iss >> bits(obj);

  return iss.good();
}

} // end of namespace wcs

#endif // __STATE_IO_IMPL_HPP__
