#include <type_traits>

namespace wcs {

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
  if constexpr (std::is_same<generator_type, std::mt19937>::value ) {
    return generator_type::state_size;
  } else if constexpr (std::is_same<generator_type, std::mt19937_64>::value ) {
    return generator_type::state_size;
  }
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

} // end of namespce wcs
