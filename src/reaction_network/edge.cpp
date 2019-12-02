#include "reaction_network/edge.hpp"

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  *  @{ */

Edge::Edge()
: m_stoichio(1), m_label("N/A")
{}

Edge::Edge(stoic_t r, const std::string& lb)
: m_stoichio(r), m_label(lb)
{}

void Edge::set_stoichiometry_ratio(stoic_t r)
{
  m_stoichio = r;
}

stoic_t Edge::get_stoichiometry_ratio() const
{
  return m_stoichio;
}

void Edge::set_label(const std::string& lb)
{
  m_label = lb;
}

std::string Edge::get_label() const
{
  return m_label;
}

/**@}*/
} // end of namespace wcs
