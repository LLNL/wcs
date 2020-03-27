#ifndef __STATE_IO_CEREAL_HPP__
#define __STATE_IO_CEREAL_HPP__

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#if defined(WCS_HAS_CEREAL)
#include <cereal/archives/binary.hpp>
#include <type_traits>
#include "state_io.hpp" // is_trivially_copyable

#define IS_CUSTOM_CEREALIZABLE(T) \
            (!std::is_arithmetic<T>::value && \
              std::is_trivially_copyable<T>::value)

namespace cereal
{
  //! Saving for trivially copyable types to binary
  template<class T> inline
  typename std::enable_if<IS_CUSTOM_CEREALIZABLE(T), void>::type
  CEREAL_SAVE_FUNCTION_NAME(BinaryOutputArchive & ar, T const & t)
  {
    ar.saveBinary(std::addressof(t), sizeof(t));
  }

  //! Loading for trivially copyable types from binary
  template<class T> inline
  typename std::enable_if<IS_CUSTOM_CEREALIZABLE(T), void>::type
  CEREAL_LOAD_FUNCTION_NAME(BinaryInputArchive & ar, T & t)
  {
    ar.loadBinary(std::addressof(t), sizeof(t));
  }
} // end of namespace cereal

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
