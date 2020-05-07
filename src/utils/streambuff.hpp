/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __STREAMBUFF_HPP__
#define __STREAMBUFF_HPP__

#include <streambuf>

namespace wcs {

/**
 * Wraps an existing buffer to use it as the internal buffer of streambuf,
 * which is then used to construct an object of basic_ostream (or one derived
 * from it). The ostream will use it as its internal streambuf.
 * Users must make sure that the external allocation outlives the object of
 * this type. The capacity of this streambuf is limited by the size of the
 * underlying buffer space, which is specified in the constructor.
 */
template<typename CharT, typename Traits = std::char_traits<CharT> >
class ostreambuff : public std::basic_streambuf<CharT, Traits>
{
 public:
  using char_type   = typename std::basic_streambuf<CharT, Traits>::char_type;
  using traits_type = typename std::basic_streambuf<CharT, Traits>::traits_type;
  using int_type    = typename std::basic_streambuf<CharT, Traits>::int_type;
  using pos_type    = typename std::basic_streambuf<CharT, Traits>::pos_type;
  using off_type    = typename std::basic_streambuf<CharT, Traits>::off_type;

  ostreambuff() = delete;

  /**
   * Set the internal buffer with the buffer of the given vector, instead of
   * relying on setbuf() or pubsetbuf(). Note that there is no other ctor
   * such that when an object of this type is created, the internal buffer
   * is set and ready.
   */
  ostreambuff(CharT* buff, size_t max_size);

  /**
   * Before the end, make sure the size of the external vector is set to the
   * exact amount of data it contains.
   */
  ~ostreambuff();

  /// Return the amount of data currently in the buffer.
  size_t size() const;

  /**
   * Return the total capacity of the underlying buffer (size of the external
   * vector).
   */
  size_t capacity() const;

  /// Shows the state of the buffer for debugging purposes
  std::ostream& print(std::ostream& os, bool show_content = false) const;

  /**
   * Limit the buffer capacity to the exact amount of data currently hold in it.
   */
  void shrink_to_fit();

 private:
  char_type* const buf;
  size_t m_capacity; ///< The maximum amount of data allowed
};


/**
 * Wraps an existing vector to use it as the internal buffer of streambuf,
 * which then is used to construct an object of basic_istream (or one derived
 * from it).
 * Users must make sure that the external vector object outlives the object of
 * this type.
 */
template<typename CharT, typename Traits = std::char_traits<CharT> >
class istreambuff : public std::basic_streambuf<CharT, Traits>
{
 public:
  using char_type   = typename std::basic_streambuf<CharT, Traits>::char_type;
  using traits_type = typename std::basic_streambuf<CharT, Traits>::traits_type;
  using int_type    = typename std::basic_streambuf<CharT, Traits>::int_type;
  using pos_type    = typename std::basic_streambuf<CharT, Traits>::pos_type;
  using off_type    = typename std::basic_streambuf<CharT, Traits>::off_type;

  istreambuff() = delete;

  istreambuff(const CharT* data, size_t sz);

  /// Return the amount of data currently in the buffer.
  size_t size() const;

  /// Shows the state of the buffer for debugging purposes
  std::ostream& print(std::ostream& os, bool show_content = false) const;


 private:
  const char_type* const buf;
  const size_t m_size; ///< The amount of data stored
};


/**
 * Wraps an existing buffer to use it as the internal buffer of streambuf,
 * which then is used to construct an object of basic_iostream (or one derived
 * from it). Users must make sure that the external buffer allocation  outlives
 * the object of this type.
 */
template<typename CharT, typename Traits = std::char_traits<CharT> >
class streambuff : public std::basic_streambuf<CharT, Traits>
{
 public:
  using char_type   = typename std::basic_streambuf<CharT, Traits>::char_type;
  using traits_type = typename std::basic_streambuf<CharT, Traits>::traits_type;
  using int_type    = typename std::basic_streambuf<CharT, Traits>::int_type;
  using pos_type    = typename std::basic_streambuf<CharT, Traits>::pos_type;
  using off_type    = typename std::basic_streambuf<CharT, Traits>::off_type;

  streambuff() = delete;

  /**
   * Set the internal buffer space with the that of the given vector, instead
   * of relying on setbuf() or pubsetbuf(). Note that there is no other ctor
   * such that when an object is created, the internal buffer is set and ready.
   * The second argument indicates whether the vector already contains data,
   * or empty.
   */
  streambuff(CharT* vec, size_t max_size, size_t cur_size = 0ul);

  /**
   * Before the end, make sure the size of the external vector is set to the
   * exact amount of data it contains.
   */
  ~streambuff();

  /// Return the amount of data currently in the buffer.
  size_t size() const;

  /**
   * Return the total capacity of the underlying buffer (size of the external
   * vector).
   */
  size_t capacity() const;

  /// Shows the state of the buffer for debugging purposes
  std::ostream& print(std::ostream& os, bool show_content = false) const;

  /**
   * Reduce the amount of memory used by the buffer to the exact amount needed.
   */
  void shrink_to_fit();

 protected:

  /**
   * Similar to the xsputn() of the base class.
   * Only writes when enough space is left in the buffer.
   */
  std::streamsize xsputn(const char_type* s, std::streamsize count) override;

  //std::streamsize xsgetn( char_type* s, std::streamsize count ) override;

 private:
  char_type* buf;
  size_t m_capacity; ///< The maximum amount of data allowed
};

} // end of namesace wcs

#include "streambuff_impl.hpp"
#endif // __STREAMBUFF_HPP__
