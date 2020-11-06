/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include <algorithm> // upper_bound
#include <cmath> // log
#include "sim_methods/ssa_sod.hpp"
#include "utils/exception.hpp"
#include "utils/seed.hpp"

#if defined(WCS_HAS_CEREAL)
#include "utils/state_io_cereal.hpp"
#endif // WCS_HAS_CEREAL

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  @{ */

SSA_SOD::SSA_SOD(const std::shared_ptr<wcs::Network>& net_ptr)
: Sim_Method(net_ptr) {}

SSA_SOD::~SSA_SOD() {}

/// Allow access to the internal random number generator for events
SSA_SOD::rng_t& SSA_SOD::rgen_e() {
  return m_rgen_evt;
}

/// Allow access to the internal random number generator for event times
SSA_SOD::rng_t& SSA_SOD::rgen_t() {
  return m_rgen_tm;
}

/**
 * Initialize the reaction propensity list by filling it with the propesity of
 * every reaction and sorting.
 */
void SSA_SOD::build_propensity_list()
{
  m_propensity.clear();

  for (const auto& vd : m_net_ptr->reaction_list())
  {
    const auto rate = m_net_ptr->get_reaction_rate(vd);
    m_propensity.emplace(priority_t{rate, static_cast<reaction_rate_t>(0.0), vd});
  }
  if (m_net_ptr->get_num_reactions() != m_propensity.size()) {
    WCS_THROW("Failed to add propensity.");
  }

  reaction_rate_t sum = static_cast<reaction_rate_t>(0.0);
  rate_idx_t::iterator it = m_propensity.begin();

  for (; it != m_propensity.end(); ++it) { // Calculate the culumative rate
    sum += it->m_rate;
  #if 0 // This does not require m_curate to be mutable
    m_propensity.modify(it, [&sum](priority_t& p) {
      p.m_curate = sum;
    });
  #else // m_curate is mutable. This avoids handling potential repositioning
    it->m_curate = sum;
  #endif
  }
}

/// Randomly determine which reaction to fire.
SSA_SOD::priority_t SSA_SOD::choose_reaction()
{
  const auto rn
    = static_cast<reaction_rate_t>(m_rgen_evt() *
                                   m_propensity.crbegin()->m_curate);

#if 0
  auto it = std::upper_bound(m_propensity.cbegin(),
                             m_propensity.cend(),
                             rn,
                             [](const double val, const priority_t& rhs)
                                -> bool { return val < rhs.m_curate; });
#else
  auto& pidx = m_propensity.get<tag_rate>();
  auto it = pidx.upper_bound(rn,
                             [](const double val, const priority_t& rhs)
                                -> bool { return val < rhs.m_curate; });
#endif
  if (it == m_propensity.cend()) {
    WCS_THROW("Failed to choose a reaction to fire");
  }
  return *it;
}

/// Randomly determine the time period until the next reaction
sim_time_t SSA_SOD::get_reaction_time()
{
  // total propensity
  const reaction_rate_t r = m_propensity.empty()?
                              static_cast<reaction_rate_t>(0) :
                              (m_propensity.crbegin()->m_curate);
  return ((r <= static_cast<reaction_rate_t>(0))?
            wcs::Network::get_etime_ulimit() :
            -static_cast<reaction_rate_t>(log(m_rgen_tm())/r));
}

/**
 * Recompute the reaction rates of those affected which are linked with
 * updating species. Also, the update cumulative propensity list.
 */
void SSA_SOD::update_reactions(const SSA_SOD::v_desc_t& vd_fired,
  const Sim_Method::affected_reactions_t& affected_reactions,
  bool check_reaction)
{
  constexpr auto zero_rate = static_cast<reaction_rate_t>(0.0);
  auto& idx_rvd = m_propensity.get<tag_rvd>();

  // update the propensity of the fired reaction
  const auto new_rate = ((check_reaction && !m_net_ptr->check_reaction(vd_fired))?
                         zero_rate : m_net_ptr->set_reaction_rate(vd_fired));
  priority_t min_new {new_rate, zero_rate, vd_fired};
  auto it_rvd = idx_rvd.find(vd_fired);
  bool ok = idx_rvd.replace(it_rvd, priority_t{new_rate, zero_rate, vd_fired});

  affected_reactions_t::const_iterator it_aff = affected_reactions.cbegin();
  // update the propensity of the rest of affected reactions
  for (; ok && (it_aff != affected_reactions.cend()); ++it_aff) {
    const auto& vd = *it_aff;
    const auto new_rate = (check_reaction && !m_net_ptr->check_reaction(vd))?
                           zero_rate : m_net_ptr->set_reaction_rate(vd);

    min_new = std::min(min_new, priority_t{new_rate, zero_rate, vd});
    auto it_rvd = idx_rvd.find(vd);
    ok = idx_rvd.replace(it_rvd, priority_t{new_rate, zero_rate, vd});
  }

  auto& idx_rate = m_propensity.get<tag_rate>();
  auto it_rate = idx_rate.find(min_new);
  if (!ok || it_rate == idx_rate.end()) {
    WCS_THROW("Failed to update reactions.");
  }

  auto it_prev = it_rate;
  reaction_rate_t sum = (it_rate != idx_rate.begin())?
                        (--it_prev)->m_curate : zero_rate;

  // update the cumulative propensity
  for (; ok && (it_rate != idx_rate.end()); ++it_rate) {
    sum += m_net_ptr->get_reaction_rate(it_rate->m_rvd);
    it_rate->m_curate = sum;
  }
  if (!ok) {
    WCS_THROW("Failed to update reactions.");
  }
}


void SSA_SOD::init(const sim_iter_t max_iter,
                   const double max_time,
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
      m_rgen_evt.set_seed();
      m_rgen_tm.set_seed();
    } else {
      seed_seq_param_t common_param_e
        = make_seed_seq_input(1, rng_seed, std::string("SSA_SOD"));
      seed_seq_param_t common_param_t
        = make_seed_seq_input(2, rng_seed, std::string("SSA_SOD"));

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

  Sim_Method::initialize_recording(m_net_ptr);

  build_propensity_list(); // prepare internal priority queue
 #if defined(WCS_HAS_ROSS)
  m_digests.emplace_back();
  m_digests.back().m_sim_time = m_sim_time;
 #endif // defined(WCS_HAS_ROSS)
}


void SSA_SOD::save_rgen_state(Sim_State_Change& digest) const
{
  constexpr size_t rng_state_size = sizeof(m_rgen_evt.engine())
                                  + sizeof(m_rgen_tm.engine());
  digest.m_rng_state.clear();
  digest.m_rng_state.reserve(rng_state_size);
  wcs::ostreamvec<char> ostrmbuf(digest.m_rng_state);
  std::ostream os(&ostrmbuf);

 #if defined(WCS_HAS_CEREAL)
  cereal::BinaryOutputArchive oarchive(os);
  oarchive(m_rgen_evt.engine(), m_rgen_tm.engine());
 #else
  os << bits(m_rgen_evt.engine()) << bits(m_rgen_tm.engine());
 #endif // defined(WCS_HAS_CEREAL)
}


void SSA_SOD::load_rgen_state(const Sim_State_Change& digest)
{
  wcs::istreamvec<char> istrmbuf(digest.m_rng_state);
  std::istream is(&istrmbuf);

 #if defined(WCS_HAS_CEREAL)
  cereal::BinaryInputArchive iarchive(is);
  iarchive(m_rgen_evt.engine(), m_rgen_tm.engine());
 #else
  is >> bits(m_rgen_evt.engine()) >> bits(m_rgen_tm.engine());
 #endif // defined(WCS_HAS_CEREAL)
}


Sim_Method::result_t SSA_SOD::schedule(sim_time_t& next_time)
{
  if (BOOST_UNLIKELY(m_propensity.empty())) { // no reaction possible
    std::cerr << "No reaction exists." << std::endl;
    return Empty;
  }

  // Determine the time when the next reaction to occur
  const auto dt = get_reaction_time();
  next_time = m_sim_time + dt;

  if (BOOST_UNLIKELY((dt >= wcs::Network::get_etime_ulimit()) ||
                     (next_time > m_max_time))) {
    std::cerr << "No more reaction can fire." << std::endl;
    return Inactive;
  }

  return Success;
}


bool SSA_SOD::forward(const sim_time_t t)
{
  if (BOOST_UNLIKELY((m_sim_iter >= m_max_iter) || (t > m_max_time))) {
    return false; // do not continue simulation
  }
  ++ m_sim_iter;
  m_sim_time = t;

 #if defined(WCS_HAS_ROSS)
  m_digests.emplace_back();
  auto& digest = m_digests.back();
  // Backup RNG state before calling choose_reaction()
  save_rgen_state(digest);
 #else
  Sim_State_Change digest;
 #endif // defined(WCS_HAS_ROSS)

  // Determine the reaction to occur at this time
  auto firing = choose_reaction();

  digest.m_sim_time = t;
  digest.m_reaction_fired = firing.m_rvd;

  // Execute the reaction, updating species counts
  Sim_Method::fire_reaction(digest);

  // Update the propensities of those reactions fired and affected
  update_reactions(digest.m_reaction_fired, digest.m_reactions_affected, true);

 #if !defined(WCS_HAS_ROSS)
  // With ROSS, tracing and sampling are moved to process at commit time
  record(firing.m_rvd);
 #endif // defined(WCS_HAS_ROSS)

  return true;
}


#if defined(WCS_HAS_ROSS)
void SSA_SOD::backward(sim_time_t& t)
{
  // State of the last event to undo
  Sim_State_Change& digest = m_digests.back();
  // The BGL vertex descriptor of the the reaction to undo
  const auto& rd_fired = digest.m_reaction_fired;

  // Undo the species update done by the reaction fired
  undo_reaction(rd_fired);
  // Undo the propensity updates done for the reactions affected
  update_reactions(rd_fired, digest.m_reactions_affected, false);

  // Restore the time
  t = digest.m_sim_time;
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


void SSA_SOD::record_first_n(const sim_iter_t num)
{
  if (m_digests.size() < 1ul) return;
  sim_iter_t i = static_cast<sim_iter_t>(0u);

  digest_list_t::iterator it = m_digests.begin();

  for (++it; it != m_digests.end(); ++it) {
    if (i >= num) break;
    record(it->m_sim_time, it->m_reaction_fired);
    i++;
  }
  m_digests.erase(m_digests.begin(), --it);
}
#endif // defined(WCS_HAS_ROSS)


std::pair<sim_iter_t, sim_time_t> SSA_SOD::run()
{
  sim_time_t t = static_cast<sim_time_t>(0);

  if (schedule(t) != Success) {
    WCS_THROW("Not able to schedule any reaction event!");
  }

  while (BOOST_LIKELY(forward(t))) {
    if (BOOST_UNLIKELY(schedule(t) != Success)) {
      break;
    }
   /*
   #if defined(WCS_HAS_ROSS)
    { // rollback test
      backward(t);
      schedule(t);
      forward(t);
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
