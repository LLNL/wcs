#ifndef __WCS_UTILS_WRITE_GRAPHVIZ_HPP__
#define __WCS_UTILS_WRITE_GRAPHVIZ_HPP__
#include <iostream>
#include <string>
#include <type_traits>


namespace wcs {
/** \addtogroup wcs_utils
 *  *  @{ */

/**
 * Check if the associative container has the interface `at()` to access the value
 * by a key. The check is done at compile time incurring no run-time overhead.
 */
template <typename Map, typename Key>
struct has_at {
private:
  typedef std::true_type yes;
  typedef std::false_type no;

  template <typename M, typename K>
  static auto test(int) ->
    decltype(std::declval<M>().at(std::declval<K>()), yes());

  template <typename M, typename K>
  static no test(...);

public:
  using type = decltype(test<Map,Key>(0));
  static constexpr bool value = std::is_same<type, yes>::value;
};


template <typename G>
std::ostream& write_graphviz(std::ostream& os, const G& g);

template <typename G>
std::ostream& operator<<(std::ostream& os, const G& g);


template <typename G>
bool write_graphviz(const std::string& out_filename, const G& g);

/**@}*/
} // end of namespace wcs

#include "write_graphviz_impl.hpp"
#endif // __WCS_UTILS_WRITE_GRAPHVIZ_HPP__
