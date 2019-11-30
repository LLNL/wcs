#include "sim_methods/sim_method.hpp"

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  *  @{ */

Sim_Method::Sim_Method()
: m_net_ptr(nullptr),
  m_max_iter(0u),
  m_enable_tracing(false),
  m_sim_time(static_cast<sim_time_t>(0)),
  m_cur_iter(0u)
{
  using directed_category
    = typename boost::graph_traits<wcs::Network::graph_t>::directed_category;

  constexpr bool is_bidirectional
    = std::is_same<directed_category, boost::bidirectional_tag>::value;

  if constexpr (!is_bidirectional) {
    WCS_THROW("Cannot get species population without in-edges.");
  }
}

Sim_Method::~Sim_Method() {}

/// Allow access to the internal tracer
Sim_Method::trace_t& Sim_Method::trace() {
  return m_trace;
}


/**@}*/
} // end of namespace wcs
