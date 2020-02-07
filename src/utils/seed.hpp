#ifndef __WCS_UTILS_SEED_HPP__
#define __WCS_UTILS_SEED_HPP__
#include <vector>
#include <unordered_set>
#include <array>
#include <random>
#include <math.h>

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

using seedseq_param_t = std::vector<uint32_t>;


/*
 *  Generate a seed sequence of N unsigned words.
 */
template <size_t N>
inline std::array<unsigned, N> compute_key(const seedseq_param_t& p)
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
                                       const seedseq_param_t& common_param,
                                       std::vector<seedseq_param_t>& unique_params)
{
  using key_t = std::template array<unsigned, N>;

  unique_params.clear();
  unique_params.reserve(num);
  if (num == 0ul) {
    return true;
  }
  else if (static_cast<double>(32*N) < log2(num)) {
    return false;
  }

  std::unordered_set<key_t> seed_set;
  uint32_t variation = 0u;
  do {
    seedseq_param_t p;
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
