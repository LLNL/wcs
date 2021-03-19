/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef WCS_SIM_METHODS_FIRE_REACTION_OMP_HPP
#define WCS_SIM_METHODS_FIRE_REACTION_OMP_HPP

#if defined(_OPENMP) && \
    defined(WCS_CACHE_DEPENDENT) && \
    defined(WCS_OMP_REACTION_REACTANTS) && \
    defined(WCS_OMP_REACTION_PRODUCTS) && \
    !defined(ENABLE_SPECIES_UPDATE_TRACKING)

#define WCS_OMP_FINE_GRAGIN 1

#include <algorithm>
#include <unordered_set>
#include "sim_methods/sim_method.hpp"


namespace wcs {
/** \addtogroup wcs_reaction_network
 *  @{ */

template <typename T>
void reduce_union(typename std::vector< std::vector<T> >& results)
{
  const size_t num_srcs = results.size();
  size_t n = results.size();
  size_t d = 1ul;

  while ((n/=2))
  {
    #pragma omp parallel for
    for (size_t i = 0ul; i < num_srcs; i += 2*d)
    {
      if (i+d >= num_srcs) {
        continue;
      }
      const typename std::vector<T>& src1 = results[i];
      const typename std::vector<T>& src2 = results[i+d];
      typename std::vector<T> uni(src1.size() + src2.size());
      typename std::vector<T>::iterator it;

      it = std::set_union(src1.begin(), src1.end(), 
                          src2.begin(), src2.end(),
                          uni.begin());
      uni.resize(it-uni.begin());
      results[i] = std::move(uni);
      results[i+d].clear();
    }
    d *= 2;
  }

  if (num_srcs % 2) {
    const typename std::vector<T>& src1 = results[0];
    const typename std::vector<T>& src2 = results[d];
    typename std::vector<T> uni(src1.size() + src2.size());
    typename std::vector<T>::iterator it;

    it = std::set_union(src1.begin(), src1.end(), 
                        src2.begin(), src2.end(),
                        uni.begin());
    uni.resize(it-uni.begin());
    results[0].swap(uni);
  }
}

/**
 * Execute the chosen reaction.
 * In other words, update the species population involved in the reaction.
 * In addition, record how the species are updated, and which other reactions
 * are affected as a result.
 * The former can be used to undo the reaction if needed. The latter is used
 * to update the propensity of the reactions affected by the changes in species
 * counts.
 */
bool Sim_Method::fire_reaction(Sim_State_Change& digest)
{
  using s_prop_t = wcs::Species;
  using v_desc_t = wcs::Network::v_desc_t;

  // Reactions do not change the connectivity, but only change the property
  // data, which are allocated outside of BGL graph, but only linked to it.
  const wcs::Network::graph_t& g = m_net_ptr->graph();
  // The vertex descriptor of the reaction to undo
  const auto& rd_firing = digest.m_reaction_fired;
  const auto& rv_firing = g[rd_firing];
  const auto& rp_firing = rv_firing.property<wcs::Network::r_prop_t>();
  const bool not_cached = rp_firing.get_dependent_reactions().empty();
  std::vector< std::vector<v_desc_t> > reactions_affected(m_num_threads);
  

 #if defined(WCS_OMP_RUN_PARTITION)
  const auto pid = m_net_ptr->get_partition_id();
 #endif // defined(WCS_OMP_RUN_PARTITION)

  // ========================= reactant species ================================
  using e_desc_t = wcs::Network::e_desc_t;
  using ie_iter_t = boost::graph_traits<Network::graph_t>::in_edge_iterator;
  ie_iter_t iei, iei_end;

  boost::tie(iei, iei_end) = boost::in_edges(rd_firing, g);
  const std::vector<e_desc_t> inedges(iei, iei_end);
  const size_t nie = inedges.size();

  #pragma omp parallel for //schedule(dynamic)
  for (size_t i = 0u; i < nie; ++i)
  {
    const auto& ei_in = inedges[i];
    const auto vd_updating = boost::source(ei_in, g);
    const auto& sv_updating = g[vd_updating];
    if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
      // in case that there are other type of vertices than species or reaction
      if (sv_updating.get_type() != wcs::Vertex::_species_) continue;
    }

    auto& sp_updating = sv_updating.property<s_prop_t>();
    const auto stoichio = g[ei_in].get_stoichiometry_ratio();
    if (stoichio == static_cast<stoic_t>(0)) {
      continue;
    }
  #ifdef NDEBUG
    sp_updating.dec_count(stoichio);
  #else
    // This really should not happen because whether the reaction is feasible is
    // checked before computing reaction time or propensity.
    if (!sp_updating.dec_count(stoichio)) { // State update
      std::string err = "Not enough reactants of " + sv_updating.get_label()
                      + "[" + std::to_string(sp_updating.get_count())
                      + "] for reaction " + g[rd_firing].get_label();
      WCS_THROW(err);
      //return false; // TODO: graceful termination
    }
  #endif

    if (not_cached)
    {
      for (const auto vi_affected :
           boost::make_iterator_range(boost::out_edges(vd_updating, g)))
      {
        const auto rd_affected = boost::target(vi_affected, g);
        if (rd_affected == rd_firing) continue;
       #if defined(WCS_OMP_RUN_PARTITION)
        if (m_net_ptr->graph()[rd_affected].get_partition() != pid) continue;
       #endif // defined(WCS_OMP_RUN_PARTITION)
        reactions_affected[omp_get_thread_num()].push_back(rd_affected);
      }
    }
  }

  // ========================== product species ================================
  using oe_iter_t = boost::graph_traits<Network::graph_t>::out_edge_iterator;
  oe_iter_t oei, oei_end;
  boost::tie(oei, oei_end) = boost::out_edges(rd_firing, g);
  const std::vector<e_desc_t> outedges(oei, oei_end);
  const size_t noe = outedges.size();

  #pragma omp parallel for //schedule(dynamic)
  for (size_t i = 0u; i < noe; ++i)
  {
    const auto& ei_out = outedges[i];
    const auto vd_updating = boost::target(ei_out, g);
    const auto& sv_updating = g[vd_updating];
    if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
      // in case that there are other type of vertices than species or reaction
      if (sv_updating.get_type() != wcs::Vertex::_species_) continue;
    }

    auto& sp_updating = sv_updating.property<s_prop_t>();
    const auto stoichio = g[ei_out].get_stoichiometry_ratio();
    if (stoichio == static_cast<stoic_t>(0)) {
      continue;
    }
  #ifdef NDEBUG
    sp_updating.inc_count(stoichio);
  #else
    if (!sp_updating.inc_count(stoichio)) { // State update
      std::string err = "Can not produce more of " + sv_updating.get_label()
                      + "[" + std::to_string(sp_updating.get_count())
                      + "] by reaction " + g[rd_firing].get_label();
      WCS_THROW(err);
      //return false; // TODO: graceful termination
    }
  #endif

    if (not_cached)
    {
      for (const auto vi_affected :
           boost::make_iterator_range(boost::out_edges(vd_updating, g)))
      {
        const auto rd_affected = boost::target(vi_affected, g);
        if (rd_affected == rd_firing) continue;
       #if defined(WCS_OMP_RUN_PARTITION)
        if (m_net_ptr->graph()[rd_affected].get_partition() != pid) continue;
       #endif // defined(_OPENMP) && defined(WCS_OMP_RUN_PARTITION)
        reactions_affected[omp_get_thread_num()].push_back(rd_affected);
      }
    }
  }

  if (not_cached) {
    reduce_union(reactions_affected);
    digest.m_dependent_reactions.assign(reactions_affected[0].begin(), reactions_affected[0].end());
  } else {
    digest.m_dependent_reactions = rp_firing.get_dependent_reactions();
  }
  return true;
}

/**@}*/
} // end of namespace wcs

#endif //defined(_OPENMP)
#endif // WCS_SIM_METHODS_FIRE_REACTION_OMP_HPP
