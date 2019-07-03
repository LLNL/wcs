#ifndef __WCS_REACTION_NETWORK_EDGE_HPP__
#define __WCS_REACTION_NETWORK_EDGE_HPP__

#include <string>
#include <iostream>
#include "wcs_types.hpp"

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  *  @{ */

class Edge {
 public:
  using stoic_t = int;

  Edge();
  Edge(stoic_t r, const std::string& lb);
  void set_stoichiometry_ratio(stoic_t r);
  stoic_t get_stoichiometry_ratio() const;
  void set_label(const std::string& lb);
  std::string get_label() const;

 protected:
  stoic_t m_stoichio; ///< stoichiometry rate
  std::string m_label; ///< label

 friend ::wcs::GraphFactory;

 template <typename G>
 friend std::ostream& ::wcs::write_graphviz_of_any_vertex_list(std::ostream& os, const G& g);
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_REACTION_NETWORK_EDGE_HPP__
