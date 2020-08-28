/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include <type_traits>
#include <cassert>
#include <algorithm>
#include "utils/exception.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

template <template <typename> typename D, typename V>
inline RNGen<D, V>::RNGen()
: m_sseq_used(false)
{
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
  if (m_sseq_used) {
    std::seed_seq sseq(m_sseq_param.begin(), m_sseq_param.end());
    m_gen.seed(sseq);
  } else {
    m_gen.seed(m_seed);
  }
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
  return m_distribution(m_gen);
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
  os << bits(m_seed) << bits(m_sseq_used) << bits(m_sseq_param) << bits(m_gen) << bits(m_distribution);
  return os;
}

template <template <typename> typename D, typename V>
template <typename S>
inline S& RNGen<D, V>::load_bits(S& is)
{
  assert (check_bits_compatibility(is));
  is >> bits(m_seed) >> bits(m_sseq_used) >> bits(m_sseq_param) >> bits(m_gen) >> bits(m_distribution);
  return is;
}

template <template <typename> typename D, typename V>
inline size_t RNGen<D, V>::byte_size() const
{
  return (sizeof(m_seed) + sizeof(m_sseq_used) +
          m_sseq_param.size() * sizeof(seed_seq_param_t::value_type) +
          sizeof(m_gen) + sizeof(m_distribution));
}

/**@}*/
} // end of namespce wcs
