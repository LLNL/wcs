/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __STREAMVEC_HPP__
#define __STREAMVEC_HPP__

#include <streambuf>

namespace wcs {

/**
 * Wraps an existing vector to use it as the internal buffer of streambuf,
 * which is then used to construct an object of basic_ostream (or one derived
 * from it). The ostream will use it as its internal streambuf.
 * Users must make sure that the external vector object outlives the object of
 * this type. This streambuf will increase the size of the underlying buffer
 * space, which is the size of the external vector, as needed.
 * In addition, the space reserving method is provided such that users can
 * preallocate the necessary space in advance to avoid the reallocation
 * overhead. Users have no way to avoid such an overhead when using stringstream
 * with a binary archive in Cereal.
 */
template<typename CharT, typename Traits = std::char_traits<CharT> >
class ostreamvec : public std::basic_streambuf<CharT, Traits>
{
 public:
  using char_type   = typename std::basic_streambuf<CharT, Traits>::char_type;
  using traits_type = typename std::basic_streambuf<CharT, Traits>::traits_type;
  using int_type    = typename std::basic_streambuf<CharT, Traits>::int_type;
  using pos_type    = typename std::basic_streambuf<CharT, Traits>::pos_type;
  using off_type    = typename std::basic_streambuf<CharT, Traits>::off_type;

  ostreamvec() = delete;

  /**
   * Set the internal buffer with the buffer of the given vector, instead of
   * relying on setbuf() or pubsetbuf(). Note that there is no other ctor
   * such that when an object of this type is created, the internal buffer
   * is set and ready.
   */
  ostreamvec(std::vector<CharT> &vec);

  /**
   * Before the end, make sure the size of the external vector is set to the
   * exact amount of data it contains.
   */
  ~ostreamvec();

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

  /**
   * Reallocate the buffer space as needed by resizing the underlying vector
   */
  void reserve(size_t n);

 protected:

  /**
   * Similar to the xsputn() of the base class except for increasing the buffer
   * capacity to accomodate whole data without failure in case of overflow.
   */
  std::streamsize xsputn( const char_type* s, std::streamsize count ) override;

  /**
   * Called inside of sputc() which is not a virtual function itself. This is to
   * make sputc() increase the internal buffer size in case of overflow.
   */
  int_type overflow(int_type c = traits_type::eof()) override;

 private:
  std::vector<char_type>& buf;
};


/**
 * Wraps an existing vector to use it as the internal buffer of streambuf,
 * which then is used to construct an object of basic_istream (or one derived
 * from it).
 * Users must make sure that the external vector object outlives the object of
 * this type.
 */
template<typename CharT, typename Traits = std::char_traits<CharT> >
class istreamvec : public std::basic_streambuf<CharT, Traits>
{
 public:
  using char_type   = typename std::basic_streambuf<CharT, Traits>::char_type;
  using traits_type = typename std::basic_streambuf<CharT, Traits>::traits_type;
  using int_type    = typename std::basic_streambuf<CharT, Traits>::int_type;
  using pos_type    = typename std::basic_streambuf<CharT, Traits>::pos_type;
  using off_type    = typename std::basic_streambuf<CharT, Traits>::off_type;

  istreamvec() = delete;

  istreamvec(const std::vector<CharT> &vec);

  /// Return the amount of data currently in the buffer.
  size_t size() const;

  /// Shows the state of the buffer for debugging purposes
  std::ostream& print(std::ostream& os, bool show_content = false) const;


 private:
  const std::vector<CharT>& buf;
};


/**
 * Wraps an existing vector to use it as the internal buffer of streambuf,
 * which then is used to construct an object of basic_iostream (or one derived
 * from it).
 * Users must make sure that the external vector object outlives the object of
 * this type. This streambuf will increase the size of the underlying buffer
 * as needed.
 */
template<typename CharT, typename Traits = std::char_traits<CharT> >
class streamvec : public std::basic_streambuf<CharT, Traits>
{
 public:
  using char_type   = typename std::basic_streambuf<CharT, Traits>::char_type;
  using traits_type = typename std::basic_streambuf<CharT, Traits>::traits_type;
  using int_type    = typename std::basic_streambuf<CharT, Traits>::int_type;
  using pos_type    = typename std::basic_streambuf<CharT, Traits>::pos_type;
  using off_type    = typename std::basic_streambuf<CharT, Traits>::off_type;

  streamvec() = delete;

  /**
   * Set the internal buffer space with the that of the given vector, instead
   * of relying on setbuf() or pubsetbuf(). Note that there is no other ctor
   * such that when an object is created, the internal buffer is set and ready.
   * The second argument indicates whether the vector already contains data,
   * or empty.
   */
  streamvec(std::vector<CharT> &vec, bool with_initial_data = false);

  /**
   * Before the end, make sure the size of the external vector is set to the
   * exact amount of data it contains.
   */
  ~streamvec();

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

  /**
   * Reallocate the buffer space as needed by resizing the underlying vector
   */
  void reserve(size_t n);

 protected:

  /**
   * Similar to the xsputn() of the base class except for increasing the buffer
   * capacity to accomodate whole data without failure in case of overflow.
   */
  std::streamsize xsputn(const char_type* s, std::streamsize count) override;

  /*
   * Similar to the xsgetn() of the base class except for updating the internal
   * pointers relevent to read operation based on the current write position.
   */
  //std::streamsize xsgetn( char_type* s, std::streamsize count ) override;

  /**
   * Called inside of sputc() which is not a virtual function itself. This is to
   * make sputc() increase the internal buffer size in case of overflow.
   */
  int_type overflow(int_type c = traits_type::eof()) override;

 private:
  std::vector<char_type>& buf;
};

} // end of namesace wcs

#include "streamvec_impl.hpp"
#endif // __STREAMVEC_HPP__
