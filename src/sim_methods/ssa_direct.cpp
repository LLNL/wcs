#include "sim_methods/ssa_direct.hpp"
#include "utils/exception.hpp"
#include <algorithm> // upper_bound
#include <cmath> // log

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  *  @{ */

SSA_Direct::SSA_Direct()
: Sim_Method() {}

SSA_Direct::~SSA_Direct() {}

/**
 * Defines the priority queue ordering by the event time of entries
 * (in ascending order).
 */
bool SSA_Direct::greater(const priority_t& v1, const priority_t& v2) {
  return (v1.first > v2.first);
}

/// Allow access to the internal random number generator
SSA_Direct::rng_t& SSA_Direct::rgen() {
  return m_rgen;
}

/**
 * Initialize the reaction propensity list by filling it with the propesity of
 * every reaction and sorting.
 */
void SSA_Direct::build_propensity_list()
{
  using r_prop_t = wcs::Reaction<v_desc_t>;

  const wcs::Network::graph_t& g = m_net_ptr->graph();

  m_propensity.clear();
  const size_t num_reactions = m_net_ptr->get_num_reactions()+1;
  m_propensity.reserve(num_reactions);
  m_pindices.reserve(num_reactions);
  constexpr auto zero_rate = static_cast<reaction_rate_t>(0.0);

  // For each reaction, check if the reaction condition is met:
  // i.e., a sufficient number of reactants considering the stoichiometry

  reaction_rate_t sum = zero_rate;
  for (const auto& vd : m_net_ptr->reaction_list())
  {
    if (m_net_ptr->check_reaction(vd)) {
      const auto& rv = g[vd]; // reaction vertex
      const auto& rp = rv.property<r_prop_t>(); // detailed vertex property data
      sum += rp.get_rate(); // cumulative reaction propensity
    }

    m_pindices.insert(std::make_pair(vd, m_propensity.size()));
    m_propensity.emplace_back(priority_t(sum, vd));
  }
  //std::stable_sort(m_propensity.begin(), m_propensity.end(), SSA_Direct::greater);
}

/// Undo the species updates applied during incomplete reaction processing.
void SSA_Direct::undo_species_updates(const std::vector<SSA_Direct::update_t>& updates) const
{
  bool ok = true;
  const wcs::Network::graph_t& g = m_net_ptr->graph();

  for (const auto& u: updates) {
    const auto& sv_undo = g[u.first];
    using s_prop_t = wcs::Species;
    auto& sp_undo = sv_undo.property<s_prop_t>();
    if (u.second > static_cast<wcs::Edge::stoic_t>(0)) {
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
    = static_cast<reaction_rate_t>(m_rgen() * m_propensity.back().first);
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

sim_time_t SSA_Direct::get_reaction_time(const SSA_Direct::priority_t& p)
{
  const reaction_rate_t r = p.first;
  return ((r <= static_cast<reaction_rate_t>(0))?
            wcs::Network::get_etime_ulimit() :
            -static_cast<reaction_rate_t>(log(m_rgen())/r));
}

/**
 * Given the reaction with the earliest time to occur, execute it (i.e.,
 * update the species population involved in the reaction). In addition,
 * record how the species are updated, and which other reactions are
 * affected as a result.
 */
bool SSA_Direct::fire_reaction(const priority_t& firing,
                            std::vector<SSA_Direct::update_t>& updating_species,
                            std::set<SSA_Direct::v_desc_t>& affected_reactions)
{
  using s_prop_t = wcs::Species;

  // Reactions do not change the connectivity, but only change the property
  // data, which are allocated outside of BGL graph, but only linked to it.
  const wcs::Network::graph_t& g = m_net_ptr->graph();

  // The BGL vertex descriptor of the curren reaction
  const auto vd_firing = firing.second;

  // reactant species
  for (const auto ei_in : boost::make_iterator_range(boost::in_edges(vd_firing, g))) {
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

    for (const auto vi_affected : boost::make_iterator_range(boost::out_edges(vd_updating, g))) {
      const auto vd_affected = boost::target(vi_affected, g);
      if (vd_affected == vd_firing) continue;
      affected_reactions.insert(vd_affected);
    }
  }

  // product species
  for (const auto ei_out : boost::make_iterator_range(boost::out_edges(vd_firing, g))) {
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

    for (const auto vi_affected : boost::make_iterator_range(boost::out_edges(vd_updating, g))) {
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
  const std::set<SSA_Direct::v_desc_t>& affected_reactions)
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
                   const unsigned rng_seed,
                   const bool enable_tracing)
{
  if (!net_ptr) {
    WCS_THROW("Invalid pointer to the reaction network.");
  }

  m_net_ptr = net_ptr;
  m_max_time = max_time;
  m_max_iter = max_iter;
  m_enable_tracing = enable_tracing;
  m_sim_time = static_cast<sim_time_t>(0);
  m_cur_iter = 0u;

  { // initialize the random number generator
    if (rng_seed == 0u)
      m_rgen.set_seed();
    else
      m_rgen.set_seed(rng_seed);

    m_rgen.param(typename rng_t::param_type(0.0, 1.0));
  }

  if (m_enable_tracing) { // record initial state of the network
    m_trace.record_initial_condition(m_net_ptr);
  }

  build_propensity_list(); // prepare internal priority queue
}

std::pair<unsigned, sim_time_t> SSA_Direct::run()
{
  // species to update as a result of the reaction fired
  std::vector<update_t> updating_species;

  // Any other reaction that takes any of species being updated as a result of
  // the current reaction as a reactant is affected. This assumes that species
  // count never reaches 'Species::m_max_count'. Otherwise, those reactions
  // that produce any of the updated species are affected as well. Here, we
  // take the assumption for simplicity. However, if we consider compartments
  // with population/concentration limits, we need to reconsider.
  std::set<v_desc_t> affected_reactions;

  for (; m_cur_iter < m_max_iter; ++ m_cur_iter) {
    if (m_propensity.empty()) { // no reaction possible
      break;
    }

    updating_species.clear();
    affected_reactions.clear();

    auto& firing = choose_reaction();
    const sim_time_t dt = get_reaction_time(firing);

    if ((dt == std::numeric_limits<sim_time_t>::infinity()) ||
        (dt >= wcs::Network::get_etime_ulimit())) {
      break;
    }

    if (!fire_reaction(firing, updating_species, affected_reactions)) {
      break;
    }

    if (m_enable_tracing) {
      m_trace.record_reaction(firing.first, firing.second);
    }

    update_reactions(firing, affected_reactions);

    m_sim_time += dt;

    if (m_sim_time >= m_max_time) {
      break;
    }
  }

  return std::make_pair(m_cur_iter, m_sim_time);
}

/**@}*/
} // end of namespace wcs
