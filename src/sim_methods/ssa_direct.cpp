#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#include <algorithm> // upper_bound
#include <cmath> // log
#include "sim_methods/ssa_direct.hpp"
#include "utils/exception.hpp"
#include "utils/seed.hpp"

#if defined(WCS_HAS_CEREAL)
#include "utils/state_io_cereal.hpp"
#endif // WCS_HAS_CEREAL

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  *  @{ */

SSA_Direct::SSA_Direct()
: Sim_Method() {}

SSA_Direct::~SSA_Direct() {}

/**
 * Defines the priority queue ordering by the event propensity
 * (in ascending order).
 */
bool SSA_Direct::less(const priority_t& v1, const priority_t& v2) {
  return (v1.first < v2.first);
}

/// Allow access to the internal random number generator for events
SSA_Direct::rng_t& SSA_Direct::rgen_e() {
  return m_rgen_evt;
}

/// Allow access to the internal random number generator for event times
SSA_Direct::rng_t& SSA_Direct::rgen_t() {
  return m_rgen_tm;
}

/**
 * Initialize the reaction propensity list by filling it with the propesity of
 * every reaction and sorting.
 */
void SSA_Direct::build_propensity_list()
{
  m_propensity.clear();
  m_pindices.clear();
  const size_t num_reactions = m_net_ptr->get_num_reactions()+1;
  m_propensity.reserve(num_reactions);
  m_pindices.reserve(num_reactions);
  constexpr auto zero_rate = static_cast<reaction_rate_t>(0.0);

  // For each reaction, check if the reaction condition is met:
  // i.e., a sufficient number of reactants considering the stoichiometry

  for (const auto& vd : m_net_ptr->reaction_list())
  {
    const auto rate = m_net_ptr->get_reaction_rate(vd);
    m_propensity.emplace_back(priority_t(rate, vd));
  }

  std::stable_sort(m_propensity.begin(), m_propensity.end(), SSA_Direct::less);

  reaction_rate_t sum = zero_rate;
  size_t i = 0ul;
  for (auto& p : m_propensity)
  {
    auto& rate = p.first;
    const auto& vd = p.second;
    sum += rate;
    rate = sum;

    m_pindices.insert(std::make_pair(vd, i++));
  }
}

/// Undo the species updates applied during incomplete reaction processing.
void SSA_Direct::undo_species_updates(
  const SSA_Direct::update_list_t& updates) const
{
  bool ok = true;
  const wcs::Network::graph_t& g = m_net_ptr->graph();

  for (const auto& u: updates) {
    const auto& sv_undo = g[u.first];
    using s_prop_t = wcs::Species;
    auto& sp_undo = sv_undo.property<s_prop_t>();
    if (u.second > static_cast<stoic_t>(0)) {
      ok &= sp_undo.dec_count(u.second);
    } else {
      ok &= sp_undo.inc_count(u.second);
    }
    if (!ok) {
      WCS_THROW("Failed to reverse the species updates");
    }
  }
}

SSA_Direct::priority_t& SSA_Direct::choose_reaction()
{
  const auto rn
    = static_cast<reaction_rate_t>(m_rgen_evt() * m_propensity.back().first);
  auto it = std::upper_bound(m_propensity.begin(),
                             m_propensity.end(),
                             rn,
                             [](const double lhs, const priority_t& rhs)
                               -> bool { return lhs < rhs.first; });
  if (it == m_propensity.end()) {
    WCS_THROW("Failed to choose a reaction to fire");
  }
  return *it;
}

sim_time_t SSA_Direct::get_reaction_time()
{
  const reaction_rate_t r = m_propensity.empty()?
                              static_cast<reaction_rate_t>(0) :
                              (m_propensity.back().first);
  return ((r <= static_cast<reaction_rate_t>(0))?
            wcs::Network::get_etime_ulimit() :
            -static_cast<reaction_rate_t>(log(m_rgen_tm())/r));
}

/**
 * Given the reaction with the earliest time to occur, execute it (i.e.,
 * update the species population involved in the reaction). In addition,
 * record how the species are updated, and which other reactions are
 * affected as a result.
 */
bool SSA_Direct::fire_reaction(const priority_t& firing,
       SSA_Direct::update_list_t& updating_species,
       SSA_Direct::affected_reactions_t& affected_reactions)
{
  using s_prop_t = wcs::Species;

  // Reactions do not change the connectivity, but only change the property
  // data, which are allocated outside of BGL graph, but only linked to it.
  const wcs::Network::graph_t& g = m_net_ptr->graph();

  // The BGL vertex descriptor of the curren reaction
  const auto vd_firing = firing.second;

  updating_species.clear();
  affected_reactions.clear();

  // reactant species
  for (const auto ei_in :
       boost::make_iterator_range(boost::in_edges(vd_firing, g)))
  {
    const auto vd_updating = boost::source(ei_in, g);
    const auto& sv_updating = g[vd_updating];
    if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
      // in case that there are other type of vertices than species or reaction
      if (sv_updating.get_type() != wcs::Vertex::_species_) continue;
    }

    auto& sp_updating = sv_updating.property<s_prop_t>();
    const auto stoichio = g[ei_in].get_stoichiometry_ratio();
    if (!sp_updating.dec_count(stoichio)) { // State update
      std::string err = "Not enough reactants of " + sv_updating.get_label()
                      + "[" + std::to_string(sp_updating.get_count())
                      + "] for reaction " + g[vd_firing].get_label();
      WCS_THROW(err);
      return false;
    }
    updating_species.emplace_back(std::make_pair(vd_updating, -stoichio));

    for (const auto vi_affected :
         boost::make_iterator_range(boost::out_edges(vd_updating, g)))
    {
      const auto vd_affected = boost::target(vi_affected, g);
      if (vd_affected == vd_firing) continue;
      affected_reactions.insert(vd_affected);
    }
  }

  // product species
  for (const auto ei_out :
       boost::make_iterator_range(boost::out_edges(vd_firing, g)))
  {
    const auto vd_updating = boost::target(ei_out, g);
    const auto& sv_updating = g[vd_updating];
    if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
      // in case that there are other type of vertices than species or reaction
      if (sv_updating.get_type() != wcs::Vertex::_species_) continue;
    }

    auto& sp_updating = sv_updating.property<s_prop_t>();
    const auto stoichio = g[ei_out].get_stoichiometry_ratio();
    if (!sp_updating.inc_count(stoichio)) { // State update
      std::string err = "Can not produce more of " + sv_updating.get_label()
                      + "[" + std::to_string(sp_updating.get_count())
                      + "] by reaction " + g[vd_firing].get_label();
      WCS_THROW(err);
      return false;
    }
    updating_species.emplace_back(std::make_pair(vd_updating, stoichio));

    for (const auto vi_affected :
         boost::make_iterator_range(boost::out_edges(vd_updating, g)))
    {
      const auto vd_affected = boost::target(vi_affected, g);
      if (vd_affected == vd_firing) continue;
      affected_reactions.insert(vd_affected);
    }
  }

  return true;
}

/**
 * Recompute the reaction rates of those affected which are linked with
 * updating species. Also, recompute the reaction time of those affected
 * and update the heap. This follows the next reaction meothod procedure.
 */
void SSA_Direct::update_reactions(priority_t& firing,
  const SSA_Direct::affected_reactions_t& affected_reactions)
{
  using r_prop_t = wcs::Reaction<v_desc_t>;

  // update the propensity of the firing reaction
  const auto vd_firing = firing.second;
  size_t pidx_min = m_pindices.at(vd_firing);
  (m_propensity.at(pidx_min)).first =  m_net_ptr->set_reaction_rate(vd_firing);

  // update the propensity of the rest of affected reactions
  for (const auto& vd : affected_reactions) {
    const size_t pidx = m_pindices.at(vd);
    pidx_min = ((pidx < pidx_min)? pidx : pidx_min);
    (m_propensity.at(pidx)).first =  m_net_ptr->set_reaction_rate(vd);
  }

  constexpr auto zero_rate = static_cast<reaction_rate_t>(0.0);
  reaction_rate_t sum = (pidx_min > 0ul)?
                        (m_propensity.at(pidx_min-1)).first : zero_rate;

  const wcs::Network::graph_t& g = m_net_ptr->graph();

  // update the cumulative propensity
  for (size_t i = pidx_min; i < m_propensity.size(); ++ i) {
    auto& prop = m_propensity.at(i);
    const auto vd = prop.second;
    const auto& rv = g[vd]; // reaction vertex
    const auto& rp = rv.property<r_prop_t>(); // detailed vertex property data
    sum += rp.get_rate(); // cumulative reaction propensity
    prop.first = sum;
  }
}


void SSA_Direct::init(std::shared_ptr<wcs::Network>& net_ptr,
                   const unsigned max_iter,
                   const double max_time,
                   const unsigned rng_seed)
{
  if (!net_ptr) {
    WCS_THROW("Invalid pointer to the reaction network.");
  }

  m_net_ptr = net_ptr;
  m_max_time = max_time;
  m_max_iter = max_iter;
  m_sim_time = static_cast<sim_time_t>(0);
  m_cur_iter = 0u;

  { // initialize the random number generator
    if (rng_seed == 0u) {
      m_rgen_evt.set_seed();
      m_rgen_tm.set_seed();
    } else {
      seed_seq_param_t common_param_e
        = make_seed_seq_input(1, rng_seed, std::string("SSA_Direct"));
      seed_seq_param_t common_param_t
        = make_seed_seq_input(2, rng_seed, std::string("SSA_Direct"));

      std::vector<seed_seq_param_t> unique_params;
      const size_t num_procs = 1ul;
      const size_t my_rank = 0ul;

      // make sure to avoid generating any duplicate seed sequence
      gen_unique_seed_seq_params<rng_t::get_state_size()>(
          num_procs, common_param_e, unique_params);
      m_rgen_evt.use_seed_seq(unique_params[my_rank]);

      // make sure to avoid generating any duplicate seed sequence
      gen_unique_seed_seq_params<rng_t::get_state_size()>(
          num_procs, common_param_t, unique_params);
      m_rgen_tm.use_seed_seq(unique_params[my_rank]);
    }

    m_rgen_evt.param(typename rng_t::param_type(0.0, 1.0));
    m_rgen_tm.param(typename rng_t::param_type(0.0, 1.0));
  }

  Sim_Method::record_initial_state(m_net_ptr);

  build_propensity_list(); // prepare internal priority queue
}

std::pair<unsigned, sim_time_t> SSA_Direct::run()
{
  // species to update as a result of the reaction fired
  update_list_t updating_species;

  // Any other reaction that takes any of species being updated as a result of
  // the current reaction as a reactant is affected. This assumes that species
  // count never reaches 'Species::m_max_count'. Otherwise, those reactions
  // that produce any of the updated species are affected as well. Here, we
  // take the assumption for simplicity. However, if we consider compartments
  // with population/concentration limits, we need to reconsider.
  affected_reactions_t affected_reactions;

  bool is_recorded = true; // initial state recording done

  for (; m_cur_iter < m_max_iter; ++ m_cur_iter) {
    if (m_propensity.empty()) { // no reaction possible
      break;
    }

    const sim_time_t dt = get_reaction_time();
    auto& firing = choose_reaction();

    if ((dt == std::numeric_limits<sim_time_t>::infinity()) ||
        (dt >= wcs::Network::get_etime_ulimit())) {
      break;
    }

    if (!fire_reaction(firing, updating_species, affected_reactions)) {
      break;
    }

    is_recorded = check_to_record(dt, firing.second);

    update_reactions(firing, affected_reactions);

    m_sim_time += dt;

    if (m_sim_time >= m_max_time) {
      if (!is_recorded) {
        record_final_state(dt, firing.second);
      }
      break;
    }
    is_recorded = false;
  }

  return std::make_pair(m_cur_iter, m_sim_time);
}

/**
 * Undo the state update done by the previous the reaction executed.
 * In addition, record which the species are restored, and which reactions are
 * affected.
 */
bool SSA_Direct::undo_reaction(const priority_t& to_undo,
       SSA_Direct::update_list_t& reverting_species,
       SSA_Direct::affected_reactions_t& affected_reactions)
{
  using s_prop_t = wcs::Species;

  // Reactions do not change the connectivity, but only change the property
  // data, which are allocated outside of BGL graph, but only linked to it.
  const wcs::Network::graph_t& g = m_net_ptr->graph();

  // The BGL vertex descriptor of the curren reaction
  const auto vd_undo = to_undo.second;

  reverting_species.clear();
  affected_reactions.clear();

  // reactant species
  for (const auto ei_in :
       boost::make_iterator_range(boost::in_edges(vd_undo, g)))
  {
    const auto vd_reverting = boost::source(ei_in, g);
    const auto& sv_reverting = g[vd_reverting];
    if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
      // in case that there are other type of vertices than species or reaction
      if (sv_reverting.get_type() != wcs::Vertex::_species_) continue;
    }

    auto& sp_reverting = sv_reverting.property<s_prop_t>();
    const auto stoichio = g[ei_in].get_stoichiometry_ratio();
    if (!sp_reverting.inc_count(stoichio)) { // State update
      std::string err = "Unable to undo the decrement of reactant "
                      + sv_reverting.get_label()
                      + "[" + std::to_string(sp_reverting.get_count())
                      + "] of reaction " + g[vd_undo].get_label();
      WCS_THROW(err);
      return false;
    }
    reverting_species.emplace_back(std::make_pair(vd_reverting, stoichio));

    for (const auto vi_affected :
         boost::make_iterator_range(boost::out_edges(vd_reverting, g)))
    {
      const auto vd_affected = boost::target(vi_affected, g);
      if (vd_affected == vd_undo) continue;
      affected_reactions.insert(vd_affected);
    }
  }

  // product species
  for (const auto ei_out :
       boost::make_iterator_range(boost::out_edges(vd_undo, g)))
  {
    const auto vd_reverting = boost::target(ei_out, g);
    const auto& sv_reverting = g[vd_reverting];
    if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
      // in case that there are other type of vertices than species or reaction
      if (sv_reverting.get_type() != wcs::Vertex::_species_) continue;
    }

    auto& sp_reverting = sv_reverting.property<s_prop_t>();
    const auto stoichio = g[ei_out].get_stoichiometry_ratio();
    if (!sp_reverting.dec_count(stoichio)) { // State update
      std::string err = "Unable to undo the production of "
                      + sv_reverting.get_label()
                      + "[" + std::to_string(sp_reverting.get_count())
                      + "] by reaction " + g[vd_undo].get_label();
      WCS_THROW(err);
      return false;
    }
    reverting_species.emplace_back(std::make_pair(vd_reverting, -stoichio));

    for (const auto vi_affected :
         boost::make_iterator_range(boost::out_edges(vd_reverting, g)))
    {
      const auto vd_affected = boost::target(vi_affected, g);
      if (vd_affected == vd_undo) continue;
      affected_reactions.insert(vd_affected);
    }
  }

  return true;
}

/**@}*/
} // end of namespace wcs
