/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

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
    auto it = this->m_idx_table.find(key); \
    if (it == this->m_idx_table.end()) { \
      WCS_THROW("The key (reaction vertex descriptor) does not exist."); \
    } \
    return it->second; \
  }; \
  \
  auto less_reaction = [this] \
    (const v_desc_t& r1, const v_desc_t& r2) -> bool \
  { \
    const auto rate1 = this->m_net_ptr->get_reaction_rate(r1); \
    const auto rate2 = this->m_net_ptr->get_reaction_rate(r2); \
    return (rate1 < rate2) || ((rate1 == rate2) && (r1 > r2)); \
  }; \
  \
  auto less_priority = [&less_reaction] \
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

SSA_NRM::SSA_NRM(const std::shared_ptr<wcs::Network>& net_ptr)
: Sim_Method(net_ptr) {}

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

sim_time_t SSA_NRM::get_reaction_time()
{
  return m_heap.front().first;
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
    return wcs::Network::get_etime_ulimit();
  }

  const auto& rv_affected = m_net_ptr->graph()[vd];
  auto& rp_affected = rv_affected.property<r_prop_t>();
  const auto rate_old = rp_affected.get_rate();
  const auto rate_new = m_net_ptr->set_reaction_rate(vd);

  if (rate_new <= static_cast<reaction_rate_t>(0)) {
    return wcs::Network::get_etime_ulimit();
  } else if (rate_old <= static_cast<reaction_rate_t>(0) ||
             rt >= wcs::Network::get_etime_ulimit()) {
    // This reaction was previously not enabled, but it is now
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

void SSA_NRM::init(const sim_iter_t max_iter,
                   const sim_time_t max_time,
                   const unsigned rng_seed)
{
  if (!m_net_ptr) {
    WCS_THROW("Invalid pointer to the reaction network.");
  }

  m_max_time = max_time;
  m_max_iter = max_iter;
  m_sim_time = static_cast<sim_time_t>(0);
  m_sim_iter = static_cast<sim_iter_t>(0u);

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

  Sim_Method::initialize_recording(m_net_ptr);

  build_heap(); // prepare internal priority queue
 #if defined(WCS_HAS_ROSS)
  m_digests.emplace_back();
  m_digests.back().m_sim_time = m_sim_time;
 #endif // defined(WCS_HAS_ROSS)
}


void SSA_NRM::save_rgen_state(Sim_State_Change& digest) const
{
  digest.m_rng_state.clear();
  digest.m_rng_state.reserve(sizeof(m_rgen.engine()));
  wcs::ostreamvec<char> ostrmbuf(digest.m_rng_state);
  std::ostream os(&ostrmbuf);

 #if defined(WCS_HAS_CEREAL)
  cereal::BinaryOutputArchive oarchive(os);
  oarchive(m_rgen.engine());
 #else
  os << bits(m_rgen_evt.engine());
 #endif // defined(WCS_HAS_CEREAL)
}


void SSA_NRM::load_rgen_state(const Sim_State_Change& digest)
{
  wcs::istreamvec<char> istrmbuf(digest.m_rng_state);
  std::istream is(&istrmbuf);

 #if defined(WCS_HAS_CEREAL)
  cereal::BinaryInputArchive iarchive(is);
  iarchive(m_rgen.engine());
 #else
  is >> bits(m_rgen.engine());
 #endif // defined(WCS_HAS_CEREAL)
}


Sim_Method::result_t SSA_NRM::schedule(revent_t& evt)
{
  if (BOOST_UNLIKELY(m_heap.empty())) { // no reaction possible
    std::cerr << "No reaction exists." << std::endl;
    return Empty;
  }

  evt = choose_reaction();

  if (BOOST_UNLIKELY(evt.first > m_max_time)) {
    std::cerr << "No more reaction can fire." << std::endl;
    return Inactive;
  }

  return Success;
}


void SSA_NRM::commit_des()
{
  if (m_recording) {
    if (m_digests.size() < 2ul) {
      WCS_THROW("No digest to commit!");
      return;
    }
    m_digests.pop_front();
    const auto& digest = m_digests.front();
    record(digest.m_sim_time, digest.m_reaction_fired);
  }
}


bool SSA_NRM::forward(const revent_t firing)
{
  const auto& t = firing.first;

  if (BOOST_UNLIKELY((m_sim_iter >= m_max_iter) || (t > m_max_time))) {
   #if defined(WCS_HAS_ROSS)
    // Whether successful or not, need to add a new digest entry such that
    // backward can remove exactly one regardless of how forward was done
    m_digests.emplace_back(firing);
   #endif // defined(WCS_HAS_ROSS)
    return false; // do not continue simulation
  }
  ++ m_sim_iter;
  m_sim_time = t;

 #if defined(WCS_HAS_ROSS)
  m_digests.emplace_back(firing);
  auto& digest = m_digests.back();

  // Backup RNG state
  save_rgen_state(digest);
  // Time of this reaction
 #else
  Sim_State_Change digest(firing);
 #endif // defined(WCS_HAS_ROSS)


  // Execute the reaction, updating species counts
  Sim_Method::fire_reaction(digest);

  // update the propensities and times of those reactions fired and affected
  update_reactions(firing, digest.m_reactions_affected, digest.m_reaction_times);

 #if !defined(WCS_HAS_ROSS)
  // With ROSS, tracing and sampling are moved to process at commit time
  record(firing.second);
 #endif // defined(WCS_HAS_ROSS)

  return true;
}


#if defined(WCS_HAS_ROSS)
void SSA_NRM::backward(revent_t& firing)
{
  // State of the last event to undo
  Sim_State_Change& digest = m_digests.back();
  // The BGL vertex descriptor of the the reaction to undo
  const auto& rd_fired = digest.m_reaction_fired;

  // Undo the species update done by the reaction fired
  undo_reaction(rd_fired);
  // Undo the propensity updates done for the reactions affected
  revert_reaction_updates(digest.m_reaction_times);

  // Restore the schedule
  firing = std::make_pair(digest.m_sim_time, digest.m_reaction_fired);
  // Restore the RNG state
  load_rgen_state(digest);
  // Free the state of the last event
  m_digests.pop_back();

  // Restore the current simulation time and iteration
  if (BOOST_UNLIKELY(m_digests.empty() ||
      (m_sim_iter == static_cast<sim_iter_t>(0)))) {
    WCS_THROW("Not able to schedule any reaction event!");
  } else {
    m_sim_time = m_digests.back().m_sim_time;
    m_sim_iter --;
  }
}


void SSA_NRM::record_first_n(const sim_iter_t num)
{
  if (m_digests.size() < 1ul) return;
  sim_iter_t i = static_cast<sim_iter_t>(0u);

  digest_list_t::iterator it = m_digests.begin();

  for (++it; it != m_digests.end(); ++it) {
    if (i >= num) break;
    record(it->m_sim_time, it->m_reaction_fired);
    i ++;
  }
  m_digests.erase(m_digests.begin(), --it);
}
#endif // defined(WCS_HAS_ROSS)


std::pair<sim_iter_t, sim_time_t> SSA_NRM::run()
{
  revent_t next_reaction;

  if (schedule(next_reaction) != Success) {
    WCS_THROW("Not able to schedule any reaction event!");
  }

  while (BOOST_LIKELY(forward(next_reaction))) {
    if (BOOST_UNLIKELY(schedule(next_reaction) != Success)) {
      break;
    }
   /*
   #if defined(WCS_HAS_ROSS)
    { // rollback test
      backward(next_reaction);
      forward(next_reaction);
      schedule(next_reaction);
    }
   #endif // defined(WCS_HAS_ROSS)
   */
  }
 #if defined(WCS_HAS_ROSS)
  record_first_n(m_sim_iter);
 #endif // defined(WCS_HAS_ROSS)

  return std::make_pair(m_sim_iter, m_sim_time);
}

/**@}*/
} // end of namespace wcs
