#ifndef TO_STRING_HPP
#define TO_STRING_HPP
#include <sstream>

namespace wcs {
template <typename T>
inline std::enable_if_t<std::is_same<T, float>::value ||
                        std::is_same<T, double>::value, std::string>
to_string_with_precision(const T a_value, const int precision = 6)
{
  std::ostringstream out;
  out.precision(precision);
  out << std::fixed << a_value;
  return out.str();
}

template <typename T>
inline std::enable_if_t<std::is_same<T, float>::value ||
                        std::is_same<T, double>::value, std::string>
to_string_in_scientific(const T a_value, const int precision = 6)
{
  std::ostringstream out;
  out.precision(precision);
  out << std::scientific << a_value;
  return out.str();
}
} // end of namespasce wcs
#endif // TO_STRING_HPP

