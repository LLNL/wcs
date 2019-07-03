#ifndef __WCS_UTILS_RNGEN_HPP__
#define __WCS_UTILS_RNGEN_HPP__
#include <random>
#include <chrono>

namespace wcs {

template <template <typename> typename D = std::uniform_real_distribution,
          typename V = double>
class RNGen {
 public:
  using result_type = V;
  using distribution_t = D<V>;
  using param_type = typename distribution_t::param_type;

  void set_seed(unsigned s);
  void set_seed();
  void param(const param_type& p);
  param_type param() const;
  result_type operator()();
  const distribution_t& distribution() const;

 protected:
  unsigned m_seed;
  std::default_random_engine m_gen;
  distribution_t m_distribution;
};

} // end of namespce wcs

#include "rngen_impl.hpp"

#endif // __WCS_UTILS_RNGEN_HPP__
