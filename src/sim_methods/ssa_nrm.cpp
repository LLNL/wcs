/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#include <algorithm>
#include "sim_methods/ssa_nrm.hpp"
#include "utils/exception.hpp"

#if defined(WCS_HAS_CEREAL)
#include "utils/state_io_cereal.hpp"
#endif // WCS_HAS_CEREAL

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  @{ */

SSA_NRM::SSA_NRM()
: Sim_Method() {}

SSA_NRM::~SSA_NRM() {}

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

/**
 * Initialize the priority queue.
 * For each reaction, compute the time to reaction using the random number
 * generator m_rgen. While doing so, check if the reaction is feasible. In case
 * that it is not, do not compute the time to reaction but set it to the
 * infinite time.
 */
void SSA_NRM::build_heap()
{
  m_heap.clear();
  m_heap.reserve(m_net_ptr->get_num_reactions()+10);
  constexpr sim_time_t unsigned_max
    = static_cast<sim_time_t>(std::numeric_limits<unsigned>::max());

  // For each reaction, check if the reaction condition is met:
  // i.e., a sufficient number of reactants

  for (const auto& vd : m_net_ptr->reaction_list())
  {
    if (!m_net_ptr->check_reaction(vd)) {
      m_heap.emplace_back(priority_t(wcs::Network::get_etime_ulimit(), vd));
    } else {
      const auto rate = m_net_ptr->get_reaction_rate(vd); // reaction rate
      const auto rn = unsigned_max/m_rgen(); // inverse of a uniform RN U(0,1)
      const auto t = log(rn)/rate;
      m_heap.emplace_back(priority_t(t, vd));
    }
  }
  std::make_heap(m_heap.begin(), m_heap.end(), SSA_NRM::later);
}

SSA_NRM::priority_t& SSA_NRM::choose_reaction()
{
  return m_heap.front();
}

sim_time_t SSA_NRM::get_reaction_time(const SSA_NRM::priority_t& p)
{
  return p.first;
}

void SSA_NRM::reset_reaction_time(const v_desc_t& vd, wcs::sim_time_t& rt)
{
  constexpr sim_time_t unsigned_max
    = static_cast<sim_time_t>(std::numeric_limits<unsigned>::max());

  const auto new_rate = m_net_ptr->set_reaction_rate(vd);

  if (!m_net_ptr->check_reaction(vd)) {
    rt = wcs::Network::get_etime_ulimit();
  } else { // Update the rate of the reaction fired
    const auto rn = unsigned_max/m_rgen();
    rt = (new_rate <= static_cast<reaction_rate_t>(0))?
          wcs::Network::get_etime_ulimit() :
          log(rn)/new_rate;

    if (std::isnan(rt)) {
      rt = wcs::Network::get_etime_ulimit();
    }
  }
}

void SSA_NRM::adjust_reaction_time(const v_desc_t& vd, wcs::sim_time_t& rt)
{
  using r_prop_t = wcs::Reaction<v_desc_t>;
  constexpr sim_time_t unsigned_max
    = static_cast<sim_time_t>(std::numeric_limits<unsigned>::max());

  const auto& rv_affected = m_net_ptr->graph()[vd];
  auto& rp_affected = rv_affected.property<r_prop_t>();
  const auto rate_old = rp_affected.get_rate();
  const auto rate_new = m_net_ptr->set_reaction_rate(vd);

  if (!m_net_ptr->check_reaction(vd)) {
    rt = wcs::Network::get_etime_ulimit();
  } else {
    if (rate_new <= static_cast<reaction_rate_t>(0)) {
      rt = wcs::Network::get_etime_ulimit();
    } else if (rate_old <= static_cast<reaction_rate_t>(0) ||
               rt >= wcs::Network::get_etime_ulimit()) {
      const auto rn = unsigned_max/m_rgen();
      rt = log(rn)/rate_new;
    } else {
      rt = rt * rate_old / rate_new;
    }

    if (std::isnan(rt)) {
      rt = wcs::Network::get_etime_ulimit();
    }
  }
}

/**
 * Recompute the reaction rates of those affected which are linked with
 * updating species. Also, recompute the reaction time of those affected
 * and update the heap. This follows the next reaction method procedure.
 */
void SSA_NRM::update_reactions(
       SSA_NRM::priority_t& firing,
       const Sim_Method::affected_reactions_t& affected,
       SSA_NRM::reaction_times_t& affected_rtimes)
{
  auto t_fired = firing.first;

  affected_rtimes.clear();
  // Store the reaction time of the the firing reaction before updating
  const std::pair<v_desc_t, sim_time_t> before_firing(firing.second, t_fired);

  reset_reaction_time(firing.second, firing.first);

  affected_reactions_t affected_reactions(affected);

  // Update the event time of the rest of affected reactions
  for (size_t hi = 1u; hi < m_heap.size(); ++ hi) {
    auto& e = m_heap[hi];
    const auto& vd = e.second; // BGL vertex descriptor of reaction node
    auto& et = e.first; // reaction time
    et -= t_fired;
    // TODO: need a better facility to locate the affected reaction entry in
    // the m_heap.
    auto it_found = affected_reactions.find(vd);
    if (it_found == affected_reactions.end()) { // not found
      continue;
    }
    // Heap items are unique. No need to find again.
    affected_reactions.erase(it_found);

    // Record the reaction time before update
    affected_rtimes.push_back(std::make_pair(vd, et));
    adjust_reaction_time(vd, et);
  }

  // Add the reaction time of the fired reaction before the update at the end
  // such that it can easily be used and removed first.
  affected_rtimes.push_back(before_firing);
  std::make_heap(m_heap.begin(), m_heap.end(), SSA_NRM::later);
}

void SSA_NRM::revert_reaction_updates(
       const wcs::sim_time_t dt,
       const SSA_NRM::reaction_times_t& affected)
{
  const auto before_firing = affected.back();

  std::map<v_desc_t, wcs::sim_time_t> affected_reactions;
  for (const auto& r : affected) {
    affected_reactions.insert(r);
  }

  // Update the event time of the rest of affected reactions
  for (size_t hi = 0u; hi < m_heap.size(); ++ hi) {
    auto& e = m_heap[hi];
    const auto& vd = e.second;
    auto& et = e.first;

    // TODO: need a better facility to locate the affected reaction entry in
    // the m_heap. red-black tree perhaps.
    auto it_found = affected_reactions.find(vd);
    if (it_found == affected_reactions.end()) { // not found
      et += dt;
    } else {
      if (vd == before_firing.first) {
        et = before_firing.second;
      } else {
        et = it_found->second + dt; // restore the previous time
      }
      m_net_ptr->set_reaction_rate(vd); // recompute reaction rate
      // Heap items are unique. No need to find again.
      affected_reactions.erase(it_found);
    }
  }

  std::make_heap(m_heap.begin(), m_heap.end(), SSA_NRM::later);
}


void SSA_NRM::init(std::shared_ptr<wcs::Network>& net_ptr,
                   const sim_iter_t max_iter,
                   const sim_time_t max_time,
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
      m_rgen.set_seed();
    } else {
      seed_seq_param_t common_param
        = make_seed_seq_input(1, rng_seed, std::string("SSA_NRM"));

      std::vector<seed_seq_param_t> unique_params;
      const size_t num_procs = 1ul;
      const size_t my_rank = 0ul;
      // make sure to avoid generating any duplicate seed sequence
      gen_unique_seed_seq_params<rng_t::get_state_size()>(
          num_procs, common_param, unique_params);
      m_rgen.use_seed_seq(unique_params[my_rank]);
    }

    constexpr unsigned uint_max = std::numeric_limits<unsigned>::max();
    m_rgen.param(typename rng_t::param_type(100, uint_max-100));
  }

  Sim_Method::record_initial_state(m_net_ptr);

  build_heap(); // prepare internal priority queue
}

std::pair<sim_iter_t, sim_time_t> SSA_NRM::run()
{
  // species to update as a result of the reaction fired
  Sim_Method::update_list_t updating_species;

  // Any other reaction that takes any of species being updated as a result of
  // the current reaction as a reactant is affected. This assumes that species
  // count never reaches 'Species::m_max_count'. Otherwise, those reactions
  // that produce any of the updated species are affected as well. Here, we
  // take the assumption for simplicity. However, if we consider compartments
  // with population/concentration limits, we need to reconsider.
  affected_reactions_t affected_reactions;
  // Backup the current reaction times of the affected reactions
  reaction_times_t reaction_times;

  bool is_recorded = true; // initial state recording done

  for (; m_cur_iter < m_max_iter; ++ m_cur_iter) {
    if (m_heap.empty()) { // no reaction possible
      break;
    }

    auto& firing = choose_reaction();
    const sim_time_t dt = get_reaction_time(firing);
    const auto vd_firing = firing.second;

    if ((dt == std::numeric_limits<sim_time_t>::infinity()) ||
        (dt >= wcs::Network::get_etime_ulimit())) {
      break;
    }

    const sim_time_t m_sim_time_last = m_sim_time;
    m_sim_time += dt;

    if (!Sim_Method::fire_reaction(vd_firing,
                                   updating_species,
                                   affected_reactions)) {
      m_sim_time = m_sim_time_last;
      break;
    }

    is_recorded = check_to_record(firing.second);

    update_reactions(firing, affected_reactions, reaction_times);

    if (m_sim_time >= m_max_time) {
      if (!is_recorded) {
        record_final_state(firing.second);
      }
      break;
    }
    is_recorded = false;
  }

  return std::make_pair(m_cur_iter, m_sim_time);
}


/**@}*/
} // end of namespace wcs
