#ifndef __WCS_UTILS_RNGEN_HPP__
#define __WCS_UTILS_RNGEN_HPP__

#include <random>
#include <chrono>

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#if defined(WCS_HAS_CEREAL)
#include <cereal/types/vector.hpp>
#include <cereal/archives/binary.hpp>
#include "utils/state_io_cereal.hpp"
ENABLE_CUSTOM_CEREAL (std::minstd_rand);
ENABLE_CUSTOM_CEREAL (std::minstd_rand0);
ENABLE_CUSTOM_CEREAL (std::mt19937)
ENABLE_CUSTOM_CEREAL (std::mt19937_64)
ENABLE_CUSTOM_CEREAL (std::uniform_int_distribution<unsigned long long>)
ENABLE_CUSTOM_CEREAL (std::uniform_int_distribution<long long>)
ENABLE_CUSTOM_CEREAL (std::uniform_int_distribution<uint32_t>)
ENABLE_CUSTOM_CEREAL (std::uniform_int_distribution<int>)
ENABLE_CUSTOM_CEREAL (std::uniform_real_distribution<double>)
ENABLE_CUSTOM_CEREAL (std::uniform_real_distribution<float>)
#endif // WCS_HAS_CEREAL

#include "utils/state_io.hpp"
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
  void use_seed_seq(const seed_seq_param_t& p);
  void param(const param_type& p);
  param_type param() const;
  result_type operator()();
  const distribution_t& distribution() const;
  //// Return the length of the generator state in words
  static constexpr unsigned get_state_size();

  /// Allow read-write acces to the internal generator engine
  generator_type& engine();
  /// Allow read-only acces to the internal generator engine
  const generator_type& engine() const;

#if defined(WCS_HAS_CEREAL)
  template <class Archive>
  void serialize( Archive & ar )
  {
    ar(m_seed, m_sseq_used, m_sseq_param, m_gen, m_distribution);
    //ar(m_seed, m_sseq_used, m_gen, m_distribution);
  }
  friend class cereal::access;
#endif // defined(WCS_HAS_CEREAL)

  template<typename S> static bool check_bits_compatibility(const S&);
  template<typename S> S& save_bits(S &os) const;
  template<typename S> S& load_bits(S &is);
  size_t byte_size() const;

 protected:
  /**
   * seed value when a single seed value is used or the master seed
   * to generate a seed sequence
   */
  unsigned m_seed;
  /// Whether to use seed_seq
  bool m_sseq_used;
  /// seed_seq input
  seed_seq_param_t m_sseq_param;
  generator_type m_gen;
  distribution_t m_distribution;
};
} // end of namespce wcs

#include "rngen_impl.hpp"

#endif // __WCS_UTILS_RNGEN_HPP__
