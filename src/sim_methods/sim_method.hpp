#ifndef __WCS_SIM_METHODS_SIM_METHOD_HPP__
#define __WCS_SIM_METHODS_SIM_METHOD_HPP__
#include <cmath>
#include <limits>
#include <unordered_map>
#include "wcs_types.hpp"
#include "reaction_network/network.hpp"
#include "utils/rngen.hpp"
#include "utils/trace_ssa.hpp"

namespace wcs {
/** \addtogroup wcs_sim_methods
 *  *  @{ */

class Sim_Method {
public:
  using v_desc_t = wcs::Network::v_desc_t;
  using trace_t = wcs::TraceSSA;
  using sim_time_t = wcs::sim_time_t;
  using reaction_rate_t = wcs::reaction_rate_t;

  Sim_Method();
  virtual ~Sim_Method();
  virtual void init(std::shared_ptr<wcs::Network>& net_ptr,
                    const unsigned max_iter,
                    const double max_time,
                    const unsigned rng_seed,
                    const bool enable_tracing) = 0;
  virtual std::pair<unsigned, sim_time_t> run() = 0;

  trace_t& trace();

protected:

  /** The pointer to the reaction network being monitored.
   *  Make sure the network object does not get destroyed
   *  while the trace refers to it.
   */
  std::shared_ptr<wcs::Network> m_net_ptr;
  unsigned m_max_iter;
  double m_max_time;
  bool m_enable_tracing;
  sim_time_t m_sim_time;
  unsigned int m_cur_iter;
  trace_t m_trace;
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_SIM_METHODS_SIM_METHOD_HPP__
