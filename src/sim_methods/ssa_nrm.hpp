#ifndef __WCS_SIM_METHODS_SSA_NRM_HPP__
#define __WCS_SIM_METHODS_SSA_NRM_HPP__
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
  using v_desc_t = wcs::Network::v_desc_t;
  using priority_t = std::pair<wcs::sim_time_t, v_desc_t>;
  using event_queue_t = std::vector<priority_t>;
  using trace_t = wcs::TraceSSA;

  /** Type for keeping track of species updates to facilitate undoing
   *  reaction processing.  */
  using update_t = std::pair<v_desc_t, wcs::Edge::stoic_t>;

  SSA_NRM();
  void init(std::shared_ptr<wcs::Network>& net_ptr,
            const unsigned max_iter,
            const double max_time,
            const unsigned rng_seed,
            const bool enable_tracing);

  std::pair<unsigned, wcs::sim_time_t> run();

  static bool later(const priority_t& v1, const priority_t& v2);
  rng_t& rgen();
  trace_t& trace();

protected:
  void build_heap();
  bool fire_reaction(priority_t& firing,
                     std::vector<update_t>& updating_species,
                     std::set<v_desc_t>& affected_reactions);
  void update_reactions(std::set<v_desc_t>& affected_reactions);
  void undo_species_updates(const std::vector<update_t>& updates) const;

protected:
  /** The pointer to the reaction network being monitored.
   *  Make sure the network object does not get destroyed
   *  while the trace refers to it.
   */
  std::shared_ptr<wcs::Network> m_net_ptr;
  unsigned m_max_iter;
  double m_max_time;
  bool m_enable_tracing;

  wcs::sim_time_t m_sim_time;
  unsigned int m_cur_iter;
  event_queue_t m_heap;
  rng_t m_rgen;
  trace_t m_trace;
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_SIM_METHODS_SSA_NRM_HPP__
