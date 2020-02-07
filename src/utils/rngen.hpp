#ifndef __WCS_UTILS_RNGEN_HPP__
#define __WCS_UTILS_RNGEN_HPP__
#include <random>
#include <chrono>
#include "utils/seed.hpp"

namespace wcs {

template <template <typename> typename D = std::uniform_real_distribution,
          typename V = double>
class RNGen {
 public:
  using result_type = V;
  using distribution_t = D<V>;
  using param_type = typename distribution_t::param_type;
  using generator_type = std::mt19937;
  /// The length of the generator state in words
  static constexpr auto state_size = generator_type::state_size;

  RNGen();

  /// Set seed when using a single value for seeding
  void set_seed(unsigned s);
  /// Set seed using the current time value
  void set_seed();
  /**
   * Set seed_seq input to generate a seed_seq object such that a sequence of
   * values (as long as the state size) rather than a single value can be used
   * for seeding
   */
  void use_seed_seq(const seedseq_param_t& p);
  void param(const param_type& p);
  param_type param() const;
  result_type operator()();
  const distribution_t& distribution() const;

 protected:
  /**
   * seed value when a single seed value is used or the master seed
   * to generate a seed sequence
   */
  unsigned m_seed;
  /// Whether to use seed_seq
  bool m_sseq_used;
  /// seed_seq input
  seedseq_param_t m_sseq_param;
  generator_type m_gen;
  distribution_t m_distribution;
};

} // end of namespce wcs

#include "rngen_impl.hpp"

#endif // __WCS_UTILS_RNGEN_HPP__
