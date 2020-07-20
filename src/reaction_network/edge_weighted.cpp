#include "reaction_network/edge_weighted.hpp"

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  @{ */

Edge_Weighted::Edge_Weighted()
: Edge(),
  m_weight(static_cast<Edge_Weighted::edge_weight_t>(0.0)),
  m_pid(unassigned_partition)
{}

Edge_Weighted::Edge_Weighted(stoic_t r, const std::string& lb)
: Edge(r, lb),
  m_weight(static_cast<Edge_Weighted::edge_weight_t>(0.0)),
  m_pid(unassigned_partition)
{}

void Edge_Weighted::set_weight(const Edge_Weighted::edge_weight_t w)
{
  m_weight = w;
}

Edge_Weighted::edge_weight_t Edge_Weighted::get_weight() const
{
  return m_weight;
}

void Edge_Weighted::set_partition(const partition_id_t pid)
{
  m_pid = pid;
}

partition_id_t Edge_Weighted::get_partition() const
{
  return m_pid;
}

/**@}*/
} // end of namespace wcs
