#ifndef WCS_TRAITS_HPP
#define WCS_TRAITS_HPP
#include <type_traits>

#if defined(__GLIBCXX__) && __GLIBCXX__ < 20150801
namespace std {
template <typename T>
struct is_trivially_copyable : integral_constant<bool, __has_trivial_copy(T)> {};
}
#endif
#endif // WCS_TRAITS_HPP
