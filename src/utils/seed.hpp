#ifndef __WCS_UTILS_SEED_HPP__
#define __WCS_UTILS_SEED_HPP__
#include <vector>
#include <unordered_set>
#include <array>
#include <random>
#include <math.h>
#include <type_traits>
#include <functional>


namespace std
{
    // hasher for std::array
    template<typename T, size_t N>
    struct hash<array<T, N> >
    {
        using arg_t = array<T, N>;
        using result_t = size_t;

        result_t operator()(const arg_t& a) const
        {
            hash<T> hasher;
            result_t h = 0ul;
            for (size_t i = 0ul; i < N; ++i)
            {
                h = h * 31 + hasher(a[i]);
            }
            return h;
        }
    };
}

namespace wcs {

/**
 * The seed_seq constructor takes an initialization list of any integer type.
 * However, we choose to use seed_seq::result_type (uint_least32_t), which is
 * the type of values `seed_seq::param()` returns.
 * The idea behind this decision is to keep the information carried in the input
 * as intact as possible. seed_seq copies the input into the internal sequence,
 * which is made of elements of at least 32 bits. How it does that depends on
 * the implementation, and may loose information if an input element relies on
 * a representation that requires more than 32 bits.
 * Users can use multiple 32-bit elements to represent such an item by using
 * `make_seed_seq_input()` provided below.
 */
using seed_seq_param_t = std::vector<std::seed_seq::result_type>;


// https://stackoverflow.com/questions/36568050/sfinae-not-happening-with-stdunderlying-type
template <typename T, bool = std::is_enum<T>::value>
struct underlying_type_SFINAE {
  using type = typename std::underlying_type<T>::type;
};

template <typename T>
struct underlying_type_SFINAE<T, false> {
  using type = T;
};


/**
 * Hash the input and return the result.
 * The internal hash algorithm we use (std::hash) generates an 64-bit integer
 * value. This interface returns a vector of 32-bit unsigned integer based on
 * it. This is to better interface with `std::seed_seq`, which internally
 * converts each input element into a value of uint_least32_t type.
 * This method can take a value of any type that `std::hash` can. Additionally,
 * an enum type is handled as `int`.
 * http://www.cplusplus.com/reference/functional/hash/
 */
template< typename T>
seed_seq_param_t make_seed_seq_input(const T& v)
{
  using item_type = std::seed_seq::result_type; // at least 32 bit

  using U = typename std::conditional<std::is_enum<T>::value,
                                      typename underlying_type_SFINAE<T>::type,
                                      T>::type;
  typename std::hash<U> h;

  const auto hv = h(static_cast<const U>(v)); // 64-bit type (i.e. size_t)

  // NOTE: This conditional is only needed if seed_seq trims a 64-bit integer
  // input element into a 32-bit one loosing the information before processing
  // it, which depends on the implementation.
  if constexpr ((sizeof(hv) == 8ul)
                && (sizeof(std::seed_seq::result_type) == 4ul)) {
    return {static_cast<item_type>(hv >> 32u),
            static_cast<item_type>(hv & 0x00000000ffffffff)};
  }
  return {static_cast<item_type>(hv)};
}

/**
 * Takes a list of arguments of heterogeneous types, each of which is not
 * necessarily an integral type, and produces a vector of 32-bit unsigned
 * integers such that a seed_seq object can be initialized with it, which
 * expects a list of integral values of a homogeneous type.
 * This method can take arguments of any type that `std::hash` can.
 * Additionally, an enum type is handled as `int`.
 */
template <typename T, typename... Args>
seed_seq_param_t make_seed_seq_input(T first, const Args&... args)
{
  seed_seq_param_t res = make_seed_seq_input(first);
  seed_seq_param_t res2 = make_seed_seq_input(args...);
  res.insert(res.end(), res2.begin(), res2.end());
  return res;
}


/*
 *  Generate a seed sequence of N unsigned words.
 */
template <size_t N>
inline std::array<unsigned, N> compute_key(const seed_seq_param_t& p)
{
  std::seed_seq ss(p.begin(), p.end());
  typename std::array<unsigned, N> gen;
  ss.generate(gen.begin(), gen.end());
  return gen;
}

/**
 *  Generates a number of seed_seq inputs, each of which consists of a
 *  variant part that is generated and the common part that is given,
 *  while making sure that the resultant seed_seq objects will be unique.
 *  The template parameter N is the state size of each seed sequence to be
 *  generated in 32-bit words. e.g., 624 for mt19937.
 */
template <size_t N>
inline bool gen_unique_seed_seq_params(const size_t num,
                                       const seed_seq_param_t& common_param,
                                       std::vector<seed_seq_param_t>& unique_params)
{
  using key_t = std::template array<unsigned, N>;
  constexpr auto bit_len = sizeof(std::seed_seq::result_type)*4*N;

  unique_params.clear();
  unique_params.reserve(num);

  if (num == 0ul) {
    return true;
  }
  else if (static_cast<double>(bit_len) < log2(num)) {
    return false;
  }

  std::unordered_set<key_t> seed_set;
  uint32_t variation = 0u;
  do {
    seed_seq_param_t p;
    p.push_back(variation);
    p.insert(p.end(), common_param.begin(), common_param.end());
    const auto s = compute_key<N>(p);
    if (seed_set.count(s) == 0ul) {
      seed_set.insert(s);
      unique_params.emplace_back(p);
    }
    variation ++;
  } while (unique_params.size() < num);
  return true;
}

} // end of namespace wcs


#endif // __WCS_UTILS_SEED_HPP__
