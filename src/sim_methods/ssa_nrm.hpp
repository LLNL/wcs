#ifndef __WCS_SSA_NRM_HPP
#define __WCS_SSA_NRM_HPP
#include <cmath>
#include <limits>
#include "wcs_types.hpp"
#include "reaction_network/network.hpp"
#include "utils/rngen.hpp"
#include "utils/trace_ssa.hpp"

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  *  @{ */

class SSA_NRM {
public:
  using rng_t = wcs::RNGen<std::uniform_int_distribution, unsigned>;
  using priority_t = std::pair<wcs::sim_time_t, wcs::Network::v_desc_t>;
  using event_queue_t = std::vector<priority_t>;
  using trace_t = wcs::TraceSSA;

  /** Type for keeping track of updating species to facilitate undoing
   *  partially done reaction processing. For example, if a reaction has
   *  two reactants and one product, executing the reaction results in
   *  decreasing the count of each reactant and increasing that of the
   *  product. However, it may be unable to decrase the second reactant
   *  count after already decreasing the that of the first one.
   *  In which case, the reaction should not fire and undo the operation.
   *  This type record which species vertex need to be increased/decreased
   *  by how much to undo the partially processed reaction.
   */
  using update_t = std::pair<wcs::Network::v_desc_t, wcs::Edge::stoic_t>;

  SSA_NRM();
  void init(std::shared_ptr<wcs::Network>& net_ptr,
            const unsigned max_iter,
            const unsigned rng_seed,
            const bool enable_tracing);
  void build_heap();
  std::pair<priority_t, bool> fire_reaction();
  std::pair<unsigned, wcs::sim_time_t> run();

  static bool later(const priority_t& v1, const priority_t& v2);
  rng_t& rgen();
  trace_t& trace();

protected:
  void undo_species_updates(const std::vector<update_t>& undos) const;

protected:
  /** The pointer to the reaction network being monitored.
   *  Make sure the network object does not get destroyed
   *  while the trace refers to it.
   */
  std::shared_ptr<wcs::Network> m_net_ptr;
  unsigned m_max_iter;
  bool m_enable_tracing;

  wcs::sim_time_t m_sim_time;
  unsigned int m_cur_iter;
  event_queue_t m_heap;
  rng_t m_rgen;
  trace_t m_trace;
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_SSA_NRM_HPP
