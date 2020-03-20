#ifndef __STATE_IO_HPP__
#define __STATE_IO_HPP__
#include "streamvec.hpp"
#include <iostream>
#include <type_traits>

#if defined(__GLIBCXX__) && __GLIBCXX__ < 20150801
namespace std {
template <typename T>
struct is_trivially_copyable : integral_constant<bool, __has_trivial_copy(T)> {};
}
#endif

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
  T v;
};

template <typename T>
typename std::enable_if<std::is_trivially_copyable<T>::value, bits_t<T&> >::type
bits(T &v)
{
  return bits_t<T&>{v};
}

template <typename T>
typename std::enable_if<std::is_trivially_copyable<T>::value, bits_t<const T&> >::type
bits(const T& v)
{
  return bits_t<const T&>{v};
}

template<typename S, typename T>
S& operator<<(S &os, const bits_t<T&>& b)
{
  os.write(reinterpret_cast<const typename S::char_type*>(&b.v), sizeof(T));
  return os;
}

template<typename S, typename T>
S& operator>>(S& is, const bits_t<T&>& b)
{
  is.read(reinterpret_cast<typename S::char_type*>(&b.v), sizeof(T));
  return is;
}


template<typename ObjT,
         typename CharT = char,
         typename Traits = std::char_traits<CharT> >
bool save_state(const ObjT& obj, std::vector<CharT>& buffer);

template<typename ObjT,
         typename CharT = char,
         typename Traits = std::char_traits<CharT> >
bool load_state(ObjT& obj, std::vector<CharT>& buffer);

} // end of namespace wcs

#include "state_io_impl.hpp"

#endif // __STATE_IO_HPP__
