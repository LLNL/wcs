/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __WCS_UTILS_RNGEN_HPP__
#define __WCS_UTILS_RNGEN_HPP__

#include <random>
#include <chrono>
#include <memory>

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#if defined(WCS_HAS_CEREAL)
#include <cereal/types/vector.hpp>
#include <cereal/types/memory.hpp>
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
/** \addtogroup wcs_utils
 *  @{ */

template <template <typename> typename D = std::uniform_real_distribution,
          typename V = double>
class RNGen {
 public:
  using result_type = V;
  using distribution_t = D<V>;
  using param_type = typename distribution_t::param_type;
  //using generator_type = std::mt19937; // better quality but has large state
  using generator_type = std::minstd_rand;
 #if WCS_THREAD_PRIVATE_RNG
  using generator_list_t = std::vector< std::unique_ptr<generator_type> >;
 #endif // WCS_THREAD_PRIVATE_RNG

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
  /**
   * This is similar to operator() in that it returns a random value drawn from
   * the current distribution. However, the difference comes from nested
   * parallel regions. The `operator()` returns a value from the thread private
   * generator identified the id of a caller thread. When the caller is not
   * a worker thread at the inner level, but the parent thread at the outer
   * level, the thread id is no longer relevant. Thefore, we pull a value from
   * the first generator object.
   */
  result_type pull();
  const distribution_t& distribution() const;
  //// Return the length of the generator state in words
  static constexpr unsigned get_state_size();

#if WCS_THREAD_PRIVATE_RNG
  /// Allow read-write acces to the internal generator engine
  generator_list_t& engine();
  /// Allow read-only acces to the internal generator engine
  const generator_list_t& engine() const;
#else
  /// Allow read-write acces to the internal generator engine
  generator_type& engine();
  /// Allow read-only acces to the internal generator engine
  const generator_type& engine() const;
#endif // WCS_THREAD_PRIVATE_RNG

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

  template<typename S> S& save_engine_bits(S &os) const;
  template<typename S> S& load_engine_bits(S &is);
  size_t engine_byte_size() const;

 #if WCS_THREAD_PRIVATE_RNG
  /**
   * Set the number of omp threads to use. By default it is set to the value
   * returned by omp_get_max_threads(). If it has to be different, call this
   * function before calling `param()`.
   */ 
  void set_num_threads(int n) { m_num_threads = n; }
  int get_num_threads() const { return m_num_threads; }
 #endif // WCS_THREAD_PRIVATE_RNG

 protected:
  using n_threads_t = uint8_t; ///< Type for number of openmp threads
  /**
   * seed value when a single seed value is used or the master seed
   * to generate a seed sequence
   */
  unsigned m_seed;
  /// Whether to use seed_seq
  bool m_sseq_used;
  /// seed_seq input
  seed_seq_param_t m_sseq_param;
#if WCS_THREAD_PRIVATE_RNG
  /// Set of thread private generator objects identifiable by the thread id
  generator_list_t m_gen;
#else
  generator_type m_gen;
#endif // WCS_THREAD_PRIVATE_RNG
  distribution_t m_distribution;

 #if WCS_THREAD_PRIVATE_RNG
  int m_num_threads;
 #endif // WCS_THREAD_PRIVATE_RNG
};

/**@}*/
} // end of namespce wcs

#include "rngen_impl.hpp"

#endif // __WCS_UTILS_RNGEN_HPP__
