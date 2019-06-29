#include "sim_methods/ssa_nrm.hpp"
#include "utils/exception.hpp"
#include <algorithm>

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  *  @{ */

SSA_NRM::SSA_NRM()
: m_net_ptr(nullptr),
  m_max_iter(0u),
  m_enable_tracing(false),
  m_sim_time(static_cast<sim_time_t>(0)),
  m_cur_iter(0u)
{
  using directed_category
    = typename boost::graph_traits<wcs::Network::graph_t>::directed_category;

  constexpr bool is_bidirectional
    = std::is_same<directed_category, boost::bidirectional_tag>::value;

  if constexpr (!is_bidirectional) {
    WCS_THROW("Cannot get species population without in-edges.");
  }
}

/**
 * Defines the priority queue ordering by the event time of entries
 * (in ascending order).
 */
bool SSA_NRM::later(const priority_t& v1, const priority_t& v2) {
  return (v1.first >= v2.first);
}

/// Allow access to the internal random number generator
SSA_NRM::rng_t& SSA_NRM::rgen() {
  return m_rgen;
}

/// Allow access to the internal tracer
SSA_NRM::trace_t& SSA_NRM::trace() {
  return m_trace;
}

/**
 * Initialize the priority queue.
 * For each reaction, compute the time to reaction using the random number
 * generator m_rgen. While doing so, check if the reaction is feasible. In case
 * that it is not, do not compute the time to reaction but set it to the
 * infinite time.
 */
void SSA_NRM::build_heap()
{
  using r_prop_t = wcs::Reaction<wcs::Network::v_desc_t>;

  const wcs::Network::graph_t& g = m_net_ptr->graph();

  m_heap.clear();
  m_heap.reserve(m_net_ptr->get_num_reactions()+10);
  constexpr double unsigned_max = static_cast<double>(std::numeric_limits<unsigned>::max());

  // For each reaction, check if the reaction condition is met:
  // i.e., a sufficient number of reactants

  for (const auto& vd : m_net_ptr->reaction_list())
  {
    if (!m_net_ptr->check_reaction(vd)) {
std::cout << "+ set the infinite time of reaction " << g[vd].get_label()  << std::endl;
      m_heap.emplace_back(priority_t(wcs::Network::get_etime_ulimit(), vd));
    } else {
      const auto& rv = g[vd]; // reaction vertex
      const auto& rp = rv.property<r_prop_t>(); // detailed vertex property data
      const auto rate = rp.get_rate(); // reaction rate
      const auto rn = unsigned_max/m_rgen();
      const auto t = log(rn)/rate;
      m_heap.emplace_back(priority_t(t, vd));
    }
  }
  std::make_heap(m_heap.begin(), m_heap.end(), SSA_NRM::later);
}

/// Undo the species updates applied during incomplete reaction processing.
void SSA_NRM::undo_species_updates(const std::vector<SSA_NRM::update_t>& undos) const
{
  const wcs::Network::graph_t& g = m_net_ptr->graph();

  for (const auto& u: undos) {
    const auto& sv_undo = g[u.first];
    using s_prop_t = wcs::Species;
    auto& sp_undo = sv_undo.property<s_prop_t>();
    if (u.second < static_cast<wcs::Edge::stoic_t>(0)) {
      sp_undo.dec_count(u.second);
    } else {
      sp_undo.inc_count(u.second);
    }
  }
}

/**
 * Pick the reaction with the earlies time to occur, execute it, update
 * the species population, recompute the reaction rates of those affected,
 * which are linked with updating species. This follows the next reaction
 * meothod procedure.
 */
std::pair<SSA_NRM::priority_t, bool> SSA_NRM::fire_reaction()
{
  using v_desc_t = wcs::Network::v_desc_t;
  using r_prop_t = wcs::Reaction<wcs::Network::v_desc_t>;
  using s_prop_t = wcs::Species;
  constexpr double unsigned_max = static_cast<double>(std::numeric_limits<unsigned>::max());

  if (m_heap.empty()) { // no reaction possible
    return std::make_pair(priority_t(), false);
  }

  // Reactions do not change the connectivity, but only change the property
  // data, which are allocated outside of BGL graph, but only linked to it.
  const wcs::Network::graph_t& g = m_net_ptr->graph();

  const auto firing = m_heap.front();
  const v_desc_t& vd_firing = firing.second; // vertex descriptor of the firing reaction
  const wcs::sim_time_t t_firing = firing.first; // time to fire the reaction
  if ((t_firing == std::numeric_limits<wcs::sim_time_t>::infinity()) ||
      (t_firing >= wcs::Network::get_etime_ulimit())) {
    return std::make_pair(firing, false);
  }

  // species to update as a result of the reaction fired
  std::vector<update_t> updating_species;
  // Any other reaction that takes any of species being updated as a result of
  // the current reaction as a reactant is affected. This assumes that species
  // count never reaches 'Species::m_max_count'. Otherwise, those reactions
  // that produce any of the updated species are affected as well. Here, we
  // take the assumption for simplicity. However, if we consider compartments
  // with population/concentration limit, we need to reconsider how to handle.
  std::set<v_desc_t> affected_reactions;

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
      // In case of an unsuccessful reaction, undo updates partially done so far
      // for this reaction. Here we first update and undo later as necessary.
      // For successful reactions, which is most likely to be the case, this
      // saves extra checking of reaction feasibility.
      // However, for optimistic PDES, a cleaner approach would be to collect
      // updating species and update all at once when it is certain that the
      // reaction can fire.
      undo_species_updates(updating_species);
      std::string err = "Not enough reactants of " + sv_updating.get_label()
                      + "[" + std::to_string(sp_updating.get_count())
                      + "] for reaction " + g[vd_firing].get_label();
      WCS_THROW(err);
      return std::make_pair(firing, false);
    }
    updating_species.emplace_back(std::make_pair(vd_updating, stoichio));

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
      // In case of an unsuccessful reaction, undo updates partially done so far
      // for this reaction.
      undo_species_updates(updating_species);
      std::string err = "Can not produce more of " + sv_updating.get_label()
                      + "[" + std::to_string(sp_updating.get_count())
                      + "] by reaction " + g[vd_firing].get_label();
      WCS_THROW(err);
      return std::make_pair(firing, false);
    }
    updating_species.emplace_back(std::make_pair(vd_updating, -stoichio));

    for (const auto vi_affected : boost::make_iterator_range(boost::out_edges(vd_updating, g))) {
      const auto vd_affected = boost::target(vi_affected, g);
      if (vd_affected == vd_firing) continue;
      affected_reactions.insert(vd_affected);
    }
  }

  if (m_enable_tracing) {
    m_trace.record_reaction(t_firing, vd_firing);
  }

  // State updates of the reaction are completed. Prepapre for the next reaction.

  auto& next_time = m_heap.front().first;
  const auto rate_new = m_net_ptr->set_reaction_rate(vd_firing);

  if (!m_net_ptr->check_reaction(vd_firing)) {
    next_time = wcs::Network::get_etime_ulimit();
  } else { // Update the rate of the reaction fired
    const auto rn = unsigned_max/m_rgen();
    next_time = (rate_new <= static_cast<wcs::sim_time_t>(0))?
                   wcs::Network::get_etime_ulimit() :
                   log(rn)/rate_new;
    if (std::isnan(next_time)) {
      next_time = wcs::Network::get_etime_ulimit();
    }
  }

  // update the event time of the rest of affected reactions
  for (size_t hi = 1u; hi < m_heap.size(); ++ hi) {
    auto& e = m_heap[hi];
    // TODO: need a better facility to locate the affected reaction entry in
    // the m_heap. red-black tree perhaps.
    auto it_found = affected_reactions.find(e.second);
    if (it_found == affected_reactions.end()) { // not found
      continue;
    }

    // Heap items are unique. No need to find again.
    affected_reactions.erase(it_found);

    auto& new_time = e.first;
    const auto& rv_affected = g[e.second];
    auto& rp_affected = rv_affected.property<r_prop_t>();
    const auto rate_old = rp_affected.get_rate();
    const auto rate_new = m_net_ptr->set_reaction_rate(e.second);

    if (!m_net_ptr->check_reaction(e.second)) {
      new_time = wcs::Network::get_etime_ulimit();
      continue; // skip reaction time computation
    }

    if (rate_new <= static_cast<wcs::sim_time_t>(0)) {
      new_time = wcs::Network::get_etime_ulimit();
    } else if (rate_old <= static_cast<wcs::sim_time_t>(0)) {
      const auto rn = unsigned_max/m_rgen();
      new_time = log(rn)/rate_new;
    } else {
      new_time = new_time * rate_old / rate_new;
    }

    if (std::isnan(new_time)) {
      new_time = wcs::Network::get_etime_ulimit();
    }
  }

  std::make_heap(m_heap.begin(), m_heap.end(), SSA_NRM::later);
  return std::make_pair(firing, true);
}


void SSA_NRM::init(std::shared_ptr<wcs::Network>& net_ptr,
                   const unsigned max_iter,
                   const unsigned rng_seed,
                   const bool enable_tracing)
{
  if (!net_ptr) {
    WCS_THROW("Invalid pointer to the reaction network.");
  }

  m_net_ptr = net_ptr;
  m_max_iter = max_iter;
  m_enable_tracing = enable_tracing;
  m_sim_time = static_cast<sim_time_t>(0);
  m_cur_iter = 0u;

  { // initialize the random number generator
    if (rng_seed == 0u)
      m_rgen.set_seed();
    else
      m_rgen.set_seed(rng_seed);

    constexpr unsigned uint_max = std::numeric_limits<unsigned>::max();
    m_rgen.param(typename rng_t::param_type(100, uint_max-100));
  }

  if (m_enable_tracing) { // record initial state of the network
    m_trace.record_initial_condition(m_net_ptr);
  }

  build_heap(); // prepare internal priority queue
}

std::pair<unsigned, wcs::sim_time_t> SSA_NRM::run()
{
  for (; m_cur_iter < m_max_iter; ++ m_cur_iter) {
    auto p = fire_reaction();
    if (!p.second) {
      break;
    }
    m_sim_time += p.first.first;
  }
  return std::make_pair(m_cur_iter, m_sim_time);
}

/**@}*/
} // end of namespace wcs
