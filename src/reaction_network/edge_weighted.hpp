#ifndef __WCS_REACTION_NETWORK_EDGE_WEIGHTED_HPP__
#define __WCS_REACTION_NETWORK_EDGE_WEIGHTED_HPP__

#include "wcs_types.hpp"
#include "reaction_network/edge.hpp"

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  @{ */

class Edge_Weighted : public Edge {
 public:
  using edge_weight_t = reaction_rate_t;

  Edge_Weighted();
  Edge_Weighted(stoic_t r, const std::string& lb);
  void set_weight(const edge_weight_t w);
  edge_weight_t get_weight() const;
  void set_partition(const partition_id_t pid);
  partition_id_t get_partition() const;

 protected:
  edge_weight_t m_weight; ///< The weight of this edge
  partition_id_t m_pid; ///< The id of the partition to which this edge belongs
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_REACTION_NETWORK_EDGE_WEIGHTED_HPP__
