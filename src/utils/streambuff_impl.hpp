/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __STREAMBUFF_IMPL_HPP__
#define __STREAMBUFF_IMPL_HPP__

#include <iostream>

namespace wcs {

//---------------------------- ostreambuff --------------------------------

template<typename CharT, typename Traits>
ostreambuff<CharT, Traits>::ostreambuff(CharT* data, size_t max_size)
: buf(data), m_capacity(max_size)
{
  static_assert(std::is_same<CharT, char_type>::value, "Invalid char_type");
  this->setp(data, data + max_size); // set pbase and epptr
}

template<typename CharT, typename Traits>
ostreambuff<CharT, Traits>::~ostreambuff()
{}

template<typename CharT, typename Traits>
size_t ostreambuff<CharT, Traits>::size() const
{
  return static_cast<size_t>(this->pptr() - this->pbase());
}

template<typename CharT, typename Traits>
size_t ostreambuff<CharT, Traits>::capacity() const
{
  return m_capacity;
}

template<typename CharT, typename Traits>
std::ostream& ostreambuff<CharT, Traits>::print(std::ostream& os,
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
void ostreambuff<CharT, Traits>::shrink_to_fit()
{
  m_capacity = size();

  this->setp(buf, buf + m_capacity); // set pbase and epptr
  this->pbump(m_capacity); // set pptr
}


//---------------------------- istreambuff --------------------------------
template<typename CharT, typename Traits>
istreambuff<CharT, Traits>::istreambuff(const CharT* data, size_t sz)
: buf(data), m_size(sz)
{
  static_assert(std::is_same<CharT, char_type>::value, "Invalid char_type");
  auto const p = const_cast<CharT*>(data);
  this->setg(p, p, p + sz); // set eback, gptr, and egptr
}

template<typename CharT, typename Traits>
size_t istreambuff<CharT, Traits>::size() const
{
  return m_size;
}

template<typename CharT, typename Traits>
std::ostream& istreambuff<CharT, Traits>::print(std::ostream& os,
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

//---------------------------- streambuff --------------------------------

template<typename CharT, typename Traits>
streambuff<CharT, Traits>::streambuff(
  CharT* data,
  size_t max_size,
  size_t cur_size)
: buf(data), m_capacity(max_size)
{
  static_assert(std::is_same<CharT, char_type>::value, "Invalid char_type");

  auto const csz = std::min(max_size, cur_size);
  auto const d_end = data + max_size;
  auto const d_cur = data + csz;

  this->setg(data, data, d_cur); // set eback, gptr, and egptr
  this->setp(data, d_end); // set pbase and epptr
  this->pbump(csz); // set pptr
}

template<typename CharT, typename Traits>
streambuff<CharT, Traits>::~streambuff()
{}

template<typename CharT, typename Traits>
size_t streambuff<CharT, Traits>::size() const
{
  return static_cast<size_t>(this->pptr() - this->pbase());
}

template<typename CharT, typename Traits>
size_t streambuff<CharT, Traits>::capacity() const
{
  return m_capacity;
}

template<typename CharT, typename Traits>
std::ostream& streambuff<CharT, Traits>::print(std::ostream& os,
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
std::streamsize streambuff<CharT, Traits>::xsputn(
  const streambuff<CharT, Traits>::char_type* s,
  std::streamsize count)
{
  if (static_cast<size_t>(count) + size() > capacity()) {
    return static_cast<std::streamsize>(0);
  }
  this->setg(this->eback(), this->gptr(), this->pptr() + count);
  return std::basic_streambuf<CharT, Traits>::xsputn(s, count);
}

/*
template<typename CharT, typename Traits>
std::streamsize streamvec<CharT, Traits>::xsgetn(
  streambuff<CharT, Traits>::char_type* s,
  std::streamsize count)
{
  this->setg(this->eback(), this->gptr(), this->pptr());
print();
  return std::basic_streambuf<CharT, Traits>::xsgetn(s, count);
}
*/

template<typename CharT, typename Traits>
void streambuff<CharT, Traits>::shrink_to_fit()
{
  const auto sz = size();
  const auto sz_read = static_cast<size_t>(this->gptr() - this->eback());

  auto const new_end  = buf + sz;
  auto const new_read = buf + std::min(sz, sz_read);

  this->setg(buf, new_read, new_end); // set eback, gptr, and egptr
  this->setp(buf, new_end); // set pbase and epptr
  this->pbump(sz); // set pptr
}

} // end of namesace wcs
#endif // __STREAMBUFF_IMPL_HPP__
