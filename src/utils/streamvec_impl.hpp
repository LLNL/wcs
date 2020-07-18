/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __STREAMVEC_IMPL_HPP__
#define __STREAMVEC_IMPL_HPP__

#include <iostream>

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

//---------------------------- ostreamvec --------------------------------

template<typename CharT, typename Traits>
ostreamvec<CharT, Traits>::ostreamvec(std::vector<CharT> &vec)
: buf(vec)
{
  static_assert(std::is_same<CharT, char_type>::value, "Invalid char_type");
  auto const p = vec.data();
  this->setp(p, p + vec.size()); // set pbase and epptr
}

template<typename CharT, typename Traits>
ostreamvec<CharT, Traits>::~ostreamvec()
{
  buf.resize(size());
}

template<typename CharT, typename Traits>
size_t ostreamvec<CharT, Traits>::size() const
{
  return static_cast<size_t>(this->pptr() - this->pbase());
}

template<typename CharT, typename Traits>
size_t ostreamvec<CharT, Traits>::capacity() const
{
  return buf.size();
}

template<typename CharT, typename Traits>
std::ostream& ostreamvec<CharT, Traits>::print(std::ostream& os,
                                               bool show_content) const
{
  const auto sz = size();
  std::stringstream ss;
  using std::operator<<;

  ss << std::hex
     << " pbase (0x" << reinterpret_cast<unsigned long long>(this->pbase())
     << ") pptr (0x"  << reinterpret_cast<unsigned long long>(this->pptr())
     << ") epptr (0x" << reinterpret_cast<unsigned long long>(this->epptr())
     << ")";

  if (show_content) {
    ss << std::endl << "buf: [";
    for (unsigned i = 0u; i < sz; ++i) {
      ss << ' ' << static_cast<unsigned int>(static_cast<unsigned char>(buf[i]));
    }
    ss << " ]";
  }

  os << "size(" << sz << ") " << ss.str() << std::endl;
  return os;
}

template<typename CharT, typename Traits>
void ostreamvec<CharT, Traits>::shrink_to_fit()
{
  const auto sz = size();
  buf.resize(sz);
  buf.shrink_to_fit();

  auto const new_base = buf.data();
  auto const new_end  = new_base + sz;

  this->setp(new_base, new_end); // set pbase and epptr
  this->pbump(sz); // set pptr
}

template<typename CharT, typename Traits>
void ostreamvec<CharT, Traits>::reserve(size_t n)
{
  const auto sz = std::min(size(), n);
  buf.resize(n);

  auto const new_base = buf.data();
  auto const new_end  = new_base + buf.size();

  this->setp(new_base, new_end); // set pbase and epptr
  this->pbump(sz); // set pptr
}

template<typename CharT, typename Traits>
std::streamsize ostreamvec<CharT, Traits>::xsputn(
  const ostreamvec<CharT, Traits>::char_type* s,
  std::streamsize count)
{
  if (static_cast<size_t>(count) + size() > capacity()) {
    reserve(size() + static_cast<size_t>(count));
  }
  return std::basic_streambuf<CharT, Traits>::xsputn(s, count);
}

template<typename CharT, typename Traits>
typename ostreamvec<CharT, Traits>::int_type
ostreamvec<CharT, Traits>::overflow(
  ostreamvec<CharT, Traits>::int_type c)
{
  reserve(size()+1ul);
  *(this->pptr()) = c;
  this->pbump(1);
  return c;
}

//---------------------------- istreamvec --------------------------------
template<typename CharT, typename Traits>
istreamvec<CharT, Traits>::istreamvec(const std::vector<CharT> &vec)
: buf(vec)
{
  static_assert(std::is_same<CharT, char_type>::value, "Invalid char_type");
  auto const p = const_cast<CharT*>(vec.data());
  this->setg(p, p, p + vec.size()); // set eback, gptr, and egptr
}

template<typename CharT, typename Traits>
size_t istreamvec<CharT, Traits>::size() const
{
  return buf.size();
}

template<typename CharT, typename Traits>
std::ostream& istreamvec<CharT, Traits>::print(std::ostream& os,
                                               bool show_content) const
{
  const auto sz = size();
  std::stringstream ss;
  using std::operator<<;

  ss << std::hex
     << " eback (0x" << reinterpret_cast<unsigned long long>(this->eback())
     << ") gptr (0x"  << reinterpret_cast<unsigned long long>(this->gptr())
     << ") egptr (0x" << reinterpret_cast<unsigned long long>(this->egptr())
     << ")";

  if (show_content) {
    ss << std::endl << "buf: [";
    for (unsigned i = 0u; i < sz; ++i) {
      ss << ' ' << static_cast<unsigned int>(static_cast<unsigned char>(buf[i]));
    }
    ss << " ]";
  }

  os << "size(" << sz << ") " << ss.str() << std::endl;
  return os;
}

//---------------------------- streamvec --------------------------------

template<typename CharT, typename Traits>
streamvec<CharT, Traits>::streamvec(
  std::vector<CharT> &vec,
  bool with_initial_data)
: buf(vec)
{
  static_assert(std::is_same<CharT, char_type>::value, "Invalid char_type");
  auto const sz = vec.size();
  auto const p = vec.data();
  auto const p_end = p + sz;

  if (with_initial_data) {
    this->setg(p, p, p_end); // set eback, gptr, and egptr
    this->setp(p, p_end); // set pbase and epptr
    this->pbump(sz); // set pptr
  } else {
    this->setg(p, p, p); // set eback, gptr, and egptr
    this->setp(p, p_end); // set pbase and epptr
  }
}

template<typename CharT, typename Traits>
streamvec<CharT, Traits>::~streamvec()
{
  buf.resize(size());
}

template<typename CharT, typename Traits>
size_t streamvec<CharT, Traits>::size() const
{
  return static_cast<size_t>(this->pptr() - this->pbase());
}

template<typename CharT, typename Traits>
size_t streamvec<CharT, Traits>::capacity() const
{
  return buf.size();
}

template<typename CharT, typename Traits>
std::ostream& streamvec<CharT, Traits>::print(std::ostream& os,
                                              bool show_content) const
{
  const auto sz = size();
  std::stringstream ss;
  using std::operator<<;

  ss << std::hex
     << " pbase (0x" << reinterpret_cast<unsigned long long>(this->pbase())
     << ") pptr (0x"  << reinterpret_cast<unsigned long long>(this->pptr())
     << ") epptr (0x" << reinterpret_cast<unsigned long long>(this->epptr())
     << ")";

  ss << std::hex
     << " eback (0x" << reinterpret_cast<unsigned long long>(this->eback())
     << ") gptr (0x"  << reinterpret_cast<unsigned long long>(this->gptr())
     << ") egptr (0x" << reinterpret_cast<unsigned long long>(this->egptr())
     << ")";

  if (show_content) {
    ss << std::endl << "buf: [";
    for (unsigned i = 0u; i < sz; ++i) {
      ss << ' ' << static_cast<unsigned int>(static_cast<unsigned char>(buf[i]));
    }
    ss << " ]";
  }

  os << "size(" << sz << ") " << ss.str() << std::endl;
  return os;
}

template<typename CharT, typename Traits>
void streamvec<CharT, Traits>::shrink_to_fit()
{
  const auto sz = size();
  const auto sz_read = static_cast<size_t>(this->gptr() - this->eback());

  buf.resize(sz);
  buf.shrink_to_fit();

  auto const new_base = buf.data();
  auto const new_end  = new_base + sz;
  auto const new_read = new_base + std::min(sz, sz_read);

  this->setg(new_base, new_read, new_end); // set eback, gptr, and egptr
  this->setp(new_base, new_end); // set pbase and epptr
  this->pbump(sz); // set pptr
}

template<typename CharT, typename Traits>
void streamvec<CharT, Traits>::reserve(size_t n)
{
  const auto sz = std::min(size(), n);
  const auto sz_read = static_cast<size_t>(this->gptr() - this->eback());

  /* If n is less than size(), it may not need to actually resize the vector
   * buf. However, by doing so, it is clear how much space is actually needed,
   * and the vector can be shriked.
   */
  buf.resize(n);

  auto const new_base = buf.data();
  auto const new_read = new_base + std::min(sz, sz_read);

  this->setg(new_base, new_read, new_base + sz); // set eback, gptr, and egptr
  this->setp(new_base, new_base + n); // set pbase and epptr
  this->pbump(sz); // set pptr
}

template<typename CharT, typename Traits>
std::streamsize streamvec<CharT, Traits>::xsputn(
  const streamvec<CharT, Traits>::char_type* s,
  std::streamsize count)
{
  if (static_cast<size_t>(count) + size() > capacity()) {
    reserve(size() + static_cast<size_t>(count));
  }
  this->setg(this->eback(), this->gptr(), this->pptr() + count);
  return std::basic_streambuf<CharT, Traits>::xsputn(s, count);
}

/*
template<typename CharT, typename Traits>
std::streamsize streamvec<CharT, Traits>::xsgetn(
  streamvec<CharT, Traits>::char_type* s,
  std::streamsize count)
{
  this->setg(this->eback(), this->gptr(), this->pptr());
print();
  return std::basic_streambuf<CharT, Traits>::xsgetn(s, count);
}
*/

template<typename CharT, typename Traits>
typename streamvec<CharT, Traits>::int_type
streamvec<CharT, Traits>::overflow(
  streamvec<CharT, Traits>::int_type c)
{
  reserve(size()+1ul);
  *(this->pptr()) = c;
  this->pbump(1);
  this->setg(this->eback(), this->gptr(), this->pptr() + 1u);
  return c;
}

/**@}*/
} // end of namesace wcs
#endif // __STREAMVEC_IMPL_HPP__
