#ifndef __WCS_REACTION_NETWORK_EDGE_WEIGHTED_HPP__
#define __WCS_REACTION_NETWORK_EDGE_WEIGHTED_HPP__

#include "wcs_types.hpp"
#include "reaction_network/edge.hpp"

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  *  @{ */

class Edge_Weighted : public Edge{
 public:

  Edge_Weighted();
  Edge_Weighted(stoic_t r, const std::string& lb);
  void set_weight(const reaction_rate_t w);
  reaction_rate_t get_weight(void) const;

 protected:
  reaction_rate_t m_weight;
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_REACTION_NETWORK_EDGE_WEIGHTED_HPP__
