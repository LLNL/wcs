namespace wcs {

template <template <typename> typename D, typename V>
inline void RNGen<D, V>::set_seed(unsigned seed)
{
  m_seed = seed;
}

template <template <typename> typename D, typename V>
inline void RNGen<D, V>::set_seed()
{
  m_seed = std::chrono::system_clock::now().time_since_epoch().count();
}

template <template <typename> typename D, typename V>
inline void RNGen<D, V>::param(const RNGen<D, V>::param_type& p)
{
  m_gen.seed(m_seed);
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

} // end of namespce wcs
