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
#include "iheap.h"

#if defined(WCS_HAS_CEREAL)
#include "utils/state_io_cereal.hpp"
#endif // WCS_HAS_CEREAL

/**
 * Due to performance, using lambdas instead of std::function or member function.
 * For each key, the indexer returns a mutable reference to the corresponding
 * entry in the in-heap index table.
 * less_reaction compares if the reaction of the first argument is likely to be
 * less frequent than that of the second one.
 * less_priority compares if the priority of the first argument is less than
 * that of the second one/
 */
#define lambdas_for_indexed_heap \
  auto indexer = [this](const v_desc_t& key) -> heap_idx_t& { \
    auto it = m_idx_table.find(key); \
    if (it == m_idx_table.end()) { \
      WCS_THROW("The key (reaction vertex descriptor) does not exist."); \
    } \
    return it->second; \
  }; \
  \
  auto less_reaction = [this] \
    (const v_desc_t& r1, const v_desc_t& r2) -> bool \
  { \
    const auto rate1 = m_net_ptr->get_reaction_rate(r1); \
    const auto rate2 = m_net_ptr->get_reaction_rate(r2); \
    return (rate1 < rate2) || ((rate1 == rate2) && (r1 > r2)); \
  }; \
  \
  auto less_priority = [this, &less_reaction] \
       (const priority_t& p1, const priority_t& p2) -> bool \
  { \
    return (p1.first > p2.first) || \
           ((p1.first == p2.first) && less_reaction(p1.second, p2.second)); \
  };
  // The previous comparator had (p1.first >= p2.first) only.
  // When comparing the result against the past commits, This needs to be
  // consistent


namespace wcs {
/** \addtogroup wcs_reaction_network
 *  @{ */

SSA_NRM::SSA_NRM()
: Sim_Method() {}

SSA_NRM::~SSA_NRM() {}

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
  lambdas_for_indexed_heap

  // For each reaction, check if the reaction condition is met:
  // i.e., a sufficient number of reactants

  m_idx_table.clear();
  m_heap.clear();
  m_heap.reserve(m_net_ptr->get_num_reactions()+10);
  constexpr sim_time_t unsigned_max
    = static_cast<sim_time_t>(std::numeric_limits<unsigned>::max());

  const auto reaction_list = m_net_ptr->reaction_list();
  for (size_t i = 0u; i < reaction_list.size(); ++i) {
    const auto& vd = reaction_list[i];

    if (!m_net_ptr->check_reaction(vd)) {
      m_heap.emplace_back(priority_t(wcs::Network::get_etime_ulimit(), vd));
    } else {
      const auto rate = m_net_ptr->get_reaction_rate(vd); // reaction rate
      const auto rn = unsigned_max/m_rgen(); // inverse of a uniform RN U(0,1)
      const auto t = log(rn)/rate;
      m_heap.emplace_back(priority_t(t, vd));
    }
    m_idx_table[vd] = static_cast<heap_idx_t>(i); // position in the heap
  }
  iheap::make(m_heap.begin(), m_heap.end(), indexer, less_priority);
}

SSA_NRM::priority_t SSA_NRM::choose_reaction()
{
  // lambdas_for_indexed_heap
  // iheap::pop(m_heap.begin(), m_heap.end(), indexer, less_priority);
  // auto p = m_heap.back();
  // m_heap.pop_back();
  // Instead of removing it and reinserting after the update,
  // leave it in the heap so as to update in place.
  return m_heap.front();
}

sim_time_t SSA_NRM::get_reaction_time(const SSA_NRM::priority_t& p)
{
  return p.first;
}

wcs::sim_time_t SSA_NRM::recompute_reaction_time(const v_desc_t& vd)
{
  constexpr sim_time_t unsigned_max
    = static_cast<sim_time_t>(std::numeric_limits<unsigned>::max());

  wcs::sim_time_t rt = wcs::Network::get_etime_ulimit();

  // Check if the reason can occur by making sure if the reactant counts are
  // at least as large as the stoichiometry, and if the product counts can
  // increase as much as the stoichiometry value.
  if (!m_net_ptr->check_reaction(vd)) {
    m_net_ptr->set_reaction_rate(vd, 0.0);
    return wcs::Network::get_etime_ulimit();
  }

  // Update the rate of the reaction fired
  const auto new_rate = m_net_ptr->set_reaction_rate(vd);

  const auto rn = unsigned_max/m_rgen();
  rt = (new_rate <= static_cast<reaction_rate_t>(0))?
        wcs::Network::get_etime_ulimit() :
        log(rn)/new_rate;

  if (std::isnan(rt)) {
    rt = wcs::Network::get_etime_ulimit();
  }

  return rt;
}

wcs::sim_time_t SSA_NRM::adjust_reaction_time(const v_desc_t& vd,
                                              wcs::sim_time_t rt)
{
  using r_prop_t = wcs::Reaction<v_desc_t>;
  constexpr sim_time_t unsigned_max
    = static_cast<sim_time_t>(std::numeric_limits<unsigned>::max());

  if (!m_net_ptr->check_reaction(vd)) {
    m_net_ptr->set_reaction_rate(vd, 0.0);
    return wcs::Network::get_etime_ulimit();
  }

  const auto& rv_affected = m_net_ptr->graph()[vd];
  auto& rp_affected = rv_affected.property<r_prop_t>();
  const auto rate_old = rp_affected.get_rate();
  const auto rate_new = m_net_ptr->set_reaction_rate(vd);

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
  return rt;
}

/**
 * Recompute the reaction rates of those affected which are linked with
 * updating species. Also, recompute the reaction time of those affected
 * and update the heap. This follows the next reaction method procedure.
 */
void SSA_NRM::update_reactions(
       const SSA_NRM::priority_t& fired,
       const Sim_Method::affected_reactions_t& affected,
       SSA_NRM::reaction_times_t& affected_rtimes)
{
  const auto& t_fired = fired.first;
  const auto& r_fired = fired.second;
  affected_rtimes.clear();

  const auto dt_fired = recompute_reaction_time(r_fired);

  lambdas_for_indexed_heap

  iheap::update(m_heap.begin(), m_heap.end(), indexer,
                r_fired, t_fired + dt_fired, less_priority);

  affected_rtimes.emplace_back(std::make_pair(r_fired, t_fired));

  for (auto& r: affected) {
    const auto t = m_heap[indexer(r)].first; // reaction time

    const auto dt = adjust_reaction_time(r, t - t_fired);
    iheap::update(m_heap.begin(), m_heap.end(), indexer,
                  r, t_fired + dt, less_priority);

    // Record the reaction time before update
    affected_rtimes.push_back(std::make_pair(r, t));
  }
}

void SSA_NRM::revert_reaction_updates(
       const wcs::sim_time_t dt,
       const SSA_NRM::reaction_times_t& affected)
{
  lambdas_for_indexed_heap

  for (auto& r: affected) {
    // Instead of recomputing the reaction rate, it could have been resotred
    // from the state saved if it was saved.
    m_net_ptr->set_reaction_rate(r.first);
    iheap::update(m_heap.begin(), m_heap.end(), indexer,
                  r.first, r.second, less_priority);
  }
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
  Sim_Method::affected_reactions_t affected_reactions;
  // Backup the current reaction times of the affected reactions
  reaction_times_t reaction_times;

  bool is_recorded = true; // initial state recording done

  for (; m_cur_iter < m_max_iter; ++ m_cur_iter) {
    if (m_heap.empty()) { // no reaction possible
      std::cerr << "No reaction exists." << std::endl;
      break;
    }

    auto firing = choose_reaction();
    const auto vd_firing = firing.second; // reaction vertex descriptor

    if (firing.first >= wcs::Network::get_etime_ulimit()) {
      std::cerr << "No more reaction can fire." << std::endl;
      break;
    }

    if (!Sim_Method::fire_reaction(vd_firing,
                                   updating_species,
                                   affected_reactions)) {
      std::cerr << "Faile to fire a reaction." << std::endl;
      break;
    }

    m_sim_time = firing.first;
    update_reactions(firing, affected_reactions, reaction_times);

    is_recorded = check_to_record(vd_firing);

    if (m_sim_time >= m_max_time) {
      if (!is_recorded) {
        record_final_state(vd_firing);
      }
      break;
    }
    is_recorded = false;
  }

  return std::make_pair(m_cur_iter, m_sim_time);
}


/**@}*/
} // end of namespace wcs
