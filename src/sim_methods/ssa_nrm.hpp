#ifndef __WCS_SIM_METHODS_SSA_NRM_HPP__
#define __WCS_SIM_METHODS_SSA_NRM_HPP__
#include <cmath>
#include <limits>
#include "sim_methods/sim_method.hpp"

namespace wcs {
/** \addtogroup wcs_sim_methods
 *  *  @{ */

class SSA_NRM : public Sim_Method {
public:
  using rng_t = wcs::RNGen<std::uniform_int_distribution, unsigned>;
  using priority_t = std::pair<wcs::sim_time_t, v_desc_t>;
  using priority_queue_t = std::vector<priority_t>;

  /** Type for keeping track of species updates to facilitate undoing
   *  reaction processing.  */
  using update_t = std::pair<v_desc_t, stoic_t>;

  SSA_NRM();
  ~SSA_NRM() override;
  void init(std::shared_ptr<wcs::Network>& net_ptr,
            const unsigned max_iter,
            const double max_time,
            const unsigned rng_seed) override;

  std::pair<unsigned, wcs::sim_time_t> run() override;

  static bool later(const priority_t& v1, const priority_t& v2);
  rng_t& rgen();

protected:
  void build_heap();
  priority_t& choose_reaction();
  sim_time_t get_reaction_time(const priority_t& p);
  bool fire_reaction(const priority_t& firing,
                     std::vector<update_t>& updating_species,
                     std::set<v_desc_t>& affected_reactions);
  void update_reactions(priority_t& firing, const std::set<v_desc_t>& affected);
  void undo_species_updates(const std::vector<update_t>& updates) const;

protected:
  priority_queue_t m_heap;
  rng_t m_rgen;
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_SIM_METHODS_SSA_NRM_HPP__
