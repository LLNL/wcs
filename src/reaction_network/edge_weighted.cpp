#include "reaction_network/edge_weighted.hpp"

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  @{ */

Edge_Weighted::Edge_Weighted()
: Edge(), m_weight(static_cast<reaction_rate_t>(0.0))
{}

Edge_Weighted::Edge_Weighted(stoic_t r, const std::string& lb)
: Edge(r, lb), m_weight(static_cast<reaction_rate_t>(0.0))
{}

void Edge_Weighted::set_weight(const reaction_rate_t w)
{
  m_weight = w;
}

reaction_rate_t Edge_Weighted::get_weight() const
{
  return m_weight;
}

/**@}*/
} // end of namespace wcs
