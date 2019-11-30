#ifndef __WCS_SIM_METHODS_SSA_DIRECT_HPP__
#define __WCS_SIM_METHODS_SSA_DIRECT_HPP__
#include <cmath>
#include <limits>
#include <unordered_map>
#include "sim_methods/sim_method.hpp"

namespace wcs {
/** \addtogroup wcs_sim_methods
 *  *  @{ */

class SSA_Direct : public Sim_Method {
public:
  using rng_t = wcs::RNGen<std::uniform_real_distribution, double>;
  using priority_t = std::pair<wcs::reaction_rate_t, v_desc_t>;
  using propensisty_list_t = std::vector<priority_t>;

  /** Type for keeping track of species updates to facilitate undoing
   *  reaction processing.  */
  using update_t = std::pair<v_desc_t, wcs::Edge::stoic_t>;

  SSA_Direct();
  ~SSA_Direct() override;
  void init(std::shared_ptr<wcs::Network>& net_ptr,
            const unsigned max_iter,
            const double max_time,
            const unsigned rng_seed,
            const bool enable_tracing) override;

  std::pair<unsigned, wcs::sim_time_t> run() override;

  static bool greater(const priority_t& v1, const priority_t& v2);
  rng_t& rgen();

protected:
  void build_propensity_list();
  priority_t& choose_reaction();
  bool fire_reaction(const priority_t& firing,
                     std::vector<update_t>& updating_species,
                     std::set<v_desc_t>& affected_reactions);
  void update_reactions(priority_t& firing, const std::set<v_desc_t>& affected);
  void undo_species_updates(const std::vector<update_t>& updates) const;

protected:
  propensisty_list_t m_propensity;
  rng_t m_rgen;
  /// map from vertex descriptor to propensity
  std::unordered_map<v_desc_t, size_t> m_pindices;
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_SIM_METHODS_SSA_DIRECT_HPP__
