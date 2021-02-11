/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include <type_traits> // std::is_same, std::is_pod
#include <cassert>     // assert
#include <algorithm>   // std::max
#include <functional>  // std::hash
#include <limits>      // std::numeric_limits
#include "utils/exception.hpp"
#include "utils/omp_diagnostics.hpp"

#if defined(WCS_HAS_OPENMP) && !defined(_OPENMP)
#error OpenMP is not enabled.
#endif

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

template <template <typename> typename D, typename V>
inline RNGen<D, V>::RNGen()
: m_sseq_used(false)
{
 #if WCS_THREAD_PRIVATE_RNG
  m_num_threads = omp_get_max_threads();
 #endif // WCS_THREAD_PRIVATE_RNG
}

template <template <typename> typename D, typename V>
inline void RNGen<D, V>::set_seed(unsigned seed)
{
  m_sseq_used = false;
  m_seed = seed;
}

template <template <typename> typename D, typename V>
inline void RNGen<D, V>::set_seed()
{
  m_sseq_used = false;
  m_seed = std::chrono::system_clock::now().time_since_epoch().count();
}

template <template <typename> typename D, typename V>
inline void RNGen<D, V>::use_seed_seq(const wcs::seed_seq_param_t& p)
{
  m_sseq_used = true;
  m_sseq_param.clear();
  m_sseq_param.assign(p.begin(), p.end());
}

template <template <typename> typename D, typename V>
inline void RNGen<D, V>::param(const RNGen<D, V>::param_type& p)
{
 #if WCS_THREAD_PRIVATE_RNG
  m_gen.resize(m_num_threads);
  assert (m_gen.size() <=
          static_cast<size_t>(std::numeric_limits<n_threads_t>::max()));

  #if OMP_DEBUG
  std::vector<wcs::my_omp_affinity> omp_aff(m_gen.size());
  #endif // OMP_DEBUG

  omp_set_dynamic(0);

  #pragma omp parallel num_threads(m_gen.size())
  {
   #if OMP_DEBUG
    omp_aff[omp_get_thread_num()].get();
   #endif // OMP_DEBUG
    const auto tid = omp_get_thread_num();
    m_gen[tid] = std::make_unique<generator_type>();
    if (m_sseq_used) {
      wcs::seed_seq_param_t sseq_thread_param;
      sseq_thread_param.reserve(m_sseq_param.size()+1);
      sseq_thread_param = m_sseq_param;
      sseq_thread_param.push_back(tid);
      std::seed_seq sseq(sseq_thread_param.begin(), sseq_thread_param.end());
      m_gen[tid]->seed(sseq);
    } else {
    // https://www.boost.org/doc/libs/1_55_0/doc/html/hash/reference.html#boost.hash_combine
      unsigned seed = m_seed ^ (std::hash<unsigned>()(tid) + 0x9e3779b9
                                + (m_seed << 6) + (m_seed >> 2));
      m_gen[tid]->seed(seed);
    }
  }

  #if OMP_DEBUG
  for (const auto& oaff: omp_aff) {
    oaff.print();
  }
  #endif // OMP_DEBUG

 #else // WCS_THREAD_PRIVATE_RNG
  if (m_sseq_used) {
    std::seed_seq sseq(m_sseq_param.begin(), m_sseq_param.end());
    m_gen.seed(sseq);
  } else {
    m_gen.seed(m_seed);
  }
 #endif // WCS_THREAD_PRIVATE_RNG
  m_distribution.param(p);
  m_distribution.reset();
}

template <template <typename> typename D, typename V>
inline typename RNGen<D, V>::param_type RNGen<D, V>::param() const
{
  return m_distribution.param();
}

template <template <typename> typename D, typename V>
inline typename RNGen<D, V>::result_type RNGen<D, V>::operator()()
{
 #if WCS_THREAD_PRIVATE_RNG
  return m_distribution(*(m_gen[omp_get_thread_num()]));
 #else
  return m_distribution(m_gen);
 #endif // WCS_THREAD_PRIVATE_RNG
}

template <template <typename> typename D, typename V>
inline typename RNGen<D, V>::result_type RNGen<D, V>::pull()
{
 #if WCS_THREAD_PRIVATE_RNG
  return m_distribution(*(m_gen[0]));
 #else
  return m_distribution(m_gen);
 #endif // WCS_THREAD_PRIVATE_RNG
}

template <template <typename> typename D, typename V>
inline const typename RNGen<D, V>::distribution_t& RNGen<D, V>::distribution() const
{
  return m_distribution;
}

template <template <typename> typename D, typename V>
constexpr unsigned RNGen<D, V>::get_state_size()
{
  if constexpr (std::is_same<generator_type, std::mt19937>::value) {
    return std::max(static_cast<unsigned>(std::mt19937::state_size), 5000u);
  } else if constexpr (std::is_same<generator_type, std::mt19937_64>::value) {
    return std::max(static_cast<unsigned>(std::mt19937_64::state_size), 2504u);
  } else if constexpr (std::is_same<generator_type, std::minstd_rand0>::value) {
    return 2u;
  } else if constexpr (std::is_same<generator_type, std::minstd_rand>::value) {
    return 2u;
  } else if constexpr (std::is_same<generator_type, std::ranlux24_base>::value) {
    return 208u;
  } else if constexpr (std::is_same<generator_type, std::ranlux24>::value) {
    return 216u;
  } else if constexpr (std::is_same<generator_type, std::ranlux48_base>::value) {
    return 112u;
  } else if constexpr (std::is_same<generator_type, std::ranlux48>::value) {
    return 120u;
  }
  // This size is only used in determining seed_seq length, and it is not
  // unsafe to use an inaccurate number for the purpose.
  return 1u;
}

#if WCS_THREAD_PRIVATE_RNG
template <template <typename> typename D, typename V>
inline typename RNGen<D, V>::generator_list_t& RNGen<D, V>::engine()
{
  return m_gen;
}

template <template <typename> typename D, typename V>
inline const typename RNGen<D, V>::generator_list_t& RNGen<D, V>::engine() const
{
  return m_gen;
}
#else // WCS_THREAD_PRIVATE_RNG
template <template <typename> typename D, typename V>
inline typename RNGen<D, V>::generator_type& RNGen<D, V>::engine()
{
  return m_gen;
}

template <template <typename> typename D, typename V>
inline const typename RNGen<D, V>::generator_type& RNGen<D, V>::engine() const
{
  return m_gen;
}
#endif // WCS_THREAD_PRIVATE_RNG

template <template <typename> typename D, typename V>
template <typename S>
inline bool RNGen<D, V>::check_bits_compatibility(const S& s)
{
  using c_t = typename S::char_type;
  using t_t = typename S::traits_type;
  const auto bp = s.rdbuf();

  if ((dynamic_cast<const streamvec<c_t, t_t>*>(bp) != nullptr) ||
      (dynamic_cast<const istreamvec<c_t, t_t>*>(bp) != nullptr) ||
      (dynamic_cast<const ostreamvec<c_t, t_t>*>(bp) != nullptr) ||
      (dynamic_cast<const streambuff<c_t, t_t>*>(bp) != nullptr) ||
      (dynamic_cast<const istreambuff<c_t, t_t>*>(bp) != nullptr) ||
      (dynamic_cast<const ostreambuff<c_t, t_t>*>(bp) != nullptr) ||
      (dynamic_cast<const std::basic_stringbuf<c_t, t_t>*>(bp) != nullptr))
  {
    return true;
  }

  return false;
}

template <template <typename> typename D, typename V>
template <typename S>
inline S& RNGen<D, V>::save_bits(S& os) const
{
  assert (check_bits_compatibility(os));
 #if WCS_THREAD_PRIVATE_RNG
  const auto num_gens = static_cast<n_threads_t>(m_gen.size());

  os << bits(m_seed) << bits(m_sseq_used) << bits(m_sseq_param)
     << bits(m_num_threads) << bits(num_gens);

  for (const auto& g: m_gen) {
    if (!!g) os << bits(*g);
  }

  os << bits(m_distribution);
 #else
  os << bits(m_seed) << bits(m_sseq_used) << bits(m_sseq_param) << bits(m_gen)
     << bits(m_distribution);
 #endif // WCS_THREAD_PRIVATE_RNG
  return os;
}

template <template <typename> typename D, typename V>
template <typename S>
inline S& RNGen<D, V>::load_bits(S& is)
{
  assert (check_bits_compatibility(is));
 #if WCS_THREAD_PRIVATE_RNG
  n_threads_t num_gens = 0u;
  is >> bits(m_seed) >> bits(m_sseq_used) >> bits(m_sseq_param)
     >> bits(m_num_threads) >> bits(num_gens);

  m_gen.resize(num_gens);

  for (auto& g: m_gen) {
     g = std::make_unique<generator_type>();
     is >> bits(*g);
  }

  is >> bits(m_distribution);
 #else
  is >> bits(m_seed) >> bits(m_sseq_used) >> bits(m_sseq_param) >> bits(m_gen)
     >> bits(m_distribution);
 #endif // WCS_THREAD_PRIVATE_RNG
  return is;
}

template <template <typename> typename D, typename V>
inline size_t RNGen<D, V>::byte_size() const
{
  // It is also assumed that the generator_type is a POD structure
  assert (std::is_pod<generator_type>::value);
 #if WCS_THREAD_PRIVATE_RNG
  assert (m_gen.size() > 0u);
  return (sizeof(m_seed) + sizeof(m_sseq_used) +
          m_sseq_param.size() * sizeof(seed_seq_param_t::value_type) +
          sizeof(n_threads_t) + sizeof(*(m_gen[0])) * m_gen.size() +
          sizeof(m_num_threads) + sizeof(m_distribution));
 #else
  return (sizeof(m_seed) + sizeof(m_sseq_used) +
          m_sseq_param.size() * sizeof(seed_seq_param_t::value_type) +
          sizeof(m_gen) + sizeof(m_distribution));
 #endif // WCS_THREAD_PRIVATE_RNG
}

template <template <typename> typename D, typename V>
template <typename S>
inline S& RNGen<D, V>::save_engine_bits(S& os) const
{
  assert (check_bits_compatibility(os));
 #if WCS_THREAD_PRIVATE_RNG
  const auto num_gens = static_cast<n_threads_t>(m_gen.size());

  os << bits(num_gens);

  for (const auto& g: m_gen) {
    if (!!g) os << bits(*g);
  }
 #else
  os << bits(m_gen);
 #endif // WCS_THREAD_PRIVATE_RNG
  return os;
}

template <template <typename> typename D, typename V>
template <typename S>
inline S& RNGen<D, V>::load_engine_bits(S& is)
{
  assert (check_bits_compatibility(is));
 #if WCS_THREAD_PRIVATE_RNG
  n_threads_t num_gens = 0u;
  is >> bits(num_gens);

  m_gen.resize(num_gens);

  for (auto& g: m_gen) {
     g = std::make_unique<generator_type>();
     is >> bits(*g);
  }
 #else
  is >> bits(m_gen);
 #endif // WCS_THREAD_PRIVATE_RNG
  return is;
}

template <template <typename> typename D, typename V>
inline size_t RNGen<D, V>::engine_byte_size() const
{
 #if WCS_THREAD_PRIVATE_RNG
  return sizeof(*(m_gen[0])) * m_gen.size();
 #else
  return sizeof(m_gen);
 #endif // WCS_THREAD_PRIVATE_RNG
}

/**@}*/
} // end of namespce wcs
