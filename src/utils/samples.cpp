/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include "utils/samples.hpp"
#include "utils/to_string.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

Samples::Samples()
: m_start_iter(static_cast<sim_iter_t>(0u)),
  m_cur_iter(static_cast<sim_iter_t>(0u)),
  m_cur_time(static_cast<sim_time_t>(0)),
  m_sample_iter_interval(static_cast<sim_iter_t>(0u)),
  m_sample_time_interval(static_cast<sim_time_t>(0)),
  m_next_sample_iter(std::numeric_limits<sim_iter_t>::max()),
  m_next_sample_time(std::numeric_limits<sim_time_t>::infinity())
{
}

void Samples::set_time_interval(const sim_time_t t_interval,
                                const sim_time_t t_start)
{
  m_sample_time_interval = t_interval;
  if (t_start > static_cast<sim_time_t>(0)) {
    m_cur_time = t_start;
  } else if (!m_samples.empty()) {
    m_cur_time = std::get<0>(m_samples.back());
  } else {
    m_cur_time = static_cast<sim_time_t>(0);
  }
  m_next_sample_time = m_cur_time + t_interval;
}

void Samples::set_iter_interval(const sim_iter_t i_interval,
                                const sim_iter_t i_start)
{
  m_sample_iter_interval = i_interval;
  m_start_iter = m_cur_iter = i_start;
  m_next_sample_iter = m_cur_iter + i_interval;
}

void Samples::record_initial_condition(const std::shared_ptr<wcs::Network>& net_ptr)
{
  m_net_ptr = net_ptr;

  if (!m_net_ptr) {
    WCS_THROW("Invaid pointer for reaction network.");
  }
  const wcs::Network::graph_t& g = m_net_ptr->graph();

  size_t i = 0ul;
  m_initial_counts.clear();
  m_initial_counts.resize(m_net_ptr->get_num_species());

  for (const auto& vd : m_net_ptr->species_list()) {
    const auto& sv = g[vd]; // vertex (property) of the species
    // detailed vertex property data of the species
    const auto& sp = sv.property<s_prop_t>();
    //m_s_id_map[vd] = i; // done in build_index()
    m_initial_counts[i++] = sp.get_count();
  }
}

void Samples::record_reaction(const sim_time_t t, const Samples::r_desc_t r)
{
  m_r_diffs[r] ++;
  m_cur_time = t;

  if (m_cur_iter ++ >= m_next_sample_iter) {
    m_next_sample_iter += m_sample_iter_interval;
    take_sample();
  } else if (m_cur_time >= m_next_sample_time) {
    m_next_sample_time += m_sample_time_interval;
    take_sample();
  }
}

void Samples::take_sample()
{
  if (!m_net_ptr) {
    WCS_THROW("Invaid pointer for reaction network.");
  }
  const wcs::Network::graph_t& g = m_net_ptr->graph();

  r_sample_t r_sample;
  s_sample_t s_sample;

  for (const auto& r_diff: m_r_diffs) {
    r_desc_t vd_reaction = r_diff.first; // BGL vertex descriptor of the reaction
    r_cnt_t rcnt = r_diff.second;
    r_sample.emplace_back(std::make_pair(vd_reaction, rcnt));

    // product species
    for (const auto ei_out :
         boost::make_iterator_range(boost::out_edges(vd_reaction, g)))
    {
      const auto vd_product = boost::target(ei_out, g);
      if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
        // in case that there are other type of vertices than species or reaction
        if (g[vd_product].get_type() != wcs::Vertex::_species_) continue;
      }
      const auto stoichio = g[ei_out].get_stoichiometry_ratio();
      m_s_diffs[vd_product] += static_cast<s_diff_t>(stoichio * rcnt);
    }

    // reactant species
    for (const auto ei_in :
         boost::make_iterator_range(boost::in_edges(vd_reaction, g)))
    {
      const auto vd_reactant = boost::source(ei_in, g);
      if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
        // in case that there are other type of vertices than species or reaction
        if (g[vd_reactant].get_type() != wcs::Vertex::_species_) continue;
      }
      const auto stoichio = g[ei_in].get_stoichiometry_ratio();
      m_s_diffs[vd_reactant] -= static_cast<s_diff_t>(stoichio * rcnt);
    }
  }

  for (const auto& s_diff: m_s_diffs) {
    s_sample.emplace_back(std::make_pair(s_diff.first, s_diff.second));
  }
  m_r_diffs.clear();
  m_s_diffs.clear();
  m_samples.emplace_back(std::make_tuple(
                           m_cur_time,
                           std::move(s_sample),
                           std::move(r_sample)));
}

/**
 * Build the map from a vertex descriptor to an index of the
 * vector for species and reaction information respectively.
 */
void Samples::build_index_maps()
{
  if (m_s_id_map.empty()) {
    size_t idx = 0ul;
    m_s_id_map.reserve(m_net_ptr->species_list().size());

    for (const auto& sd : m_net_ptr->species_list()) {
      m_s_id_map[sd] = idx++;
    }
  }

  if (m_r_id_map.empty()) {
    size_t idx = 0ul;
    m_r_id_map.reserve(m_net_ptr->reaction_list().size());

    for (const auto& rd : m_net_ptr->reaction_list()) {
      m_r_id_map[rd] = idx++;
    }
  }
}

/**
 * Write the header (species labels), and write the initial species population.
 */
std::ostream& Samples::write_header(std::ostream& os, size_t num_reactions) const
{ // write the header to show the species labels and the initial population
  const auto num_species = m_net_ptr->get_num_species();

  sim_iter_t m_num_events = m_cur_iter - m_start_iter;

  std::string ostr = "num_species = " + std::to_string(num_species)
                   + "\tnum_reactions = " + std::to_string(num_reactions)
                   + "\tnum_events = " + std::to_string(m_num_events)
                   + "\nTime: ";

  ostr.reserve(ostr.size() + num_species*cnt_digits + num_reactions*2 + 1);

  ostr += m_net_ptr->show_species_labels("")
        + m_net_ptr->show_reaction_labels("")
        + '\n' + std::to_string(0.000000);

  // write the initial population
  for (const auto scnt : m_initial_counts) {
    ostr += '\t' + std::to_string(scnt);
  }
  // write the initial reaction distribution
  for (size_t i = 0ul; i < num_reactions; i++) {
    ostr += "\t0";
  }
  os << ostr << '\n';
  return os;
}

void Samples::count_species(
  const Samples::s_sample_t& ss,
  std::vector<species_cnt_t>& species) const
{
  for (const auto& s: ss) {
    species.at(m_s_id_map.at(std::get<0>(s))) += std::get<1>(s);
  }
}

void Samples::count_reactions(
  const Samples::r_sample_t& rs,
  std::vector<Samples::r_cnt_t>& reactions) const
{
  for (const auto& r: rs) {
    reactions.at(m_r_id_map.at(std::get<0>(r))) += std::get<1>(r);
  }
}

std::ostream& Samples::print_stats(
  const sim_time_t sim_time,
  const std::vector<species_cnt_t>& species,
  const std::vector<Samples::r_cnt_t>& reactions,
  std::string& tmpstr,
  std::ostream& os) const
{
  tmpstr = to_string_in_scientific(sim_time);
  const size_t tstr_sz = tmpstr.size();
  tmpstr.reserve(tstr_sz + 2 +
                 species.size()*cnt_digits + reactions.size()*cnt_digits);

  for (const auto scnt : species) {
    tmpstr += '\t' + std::to_string(scnt);
  }
  for (const auto rcnt : reactions) {
    tmpstr += '\t' + std::to_string(rcnt);
  }
  os << tmpstr << '\n';
  return os;
}

std::ostream& Samples::write(std::ostream& os)
{
  // species population state
  std::vector<species_cnt_t> species(m_initial_counts);

  /// Show how many times each reaction fires
  std::vector<r_cnt_t> reactions;
  reactions.resize(m_net_ptr->reaction_list().size(), 0ul);

  std::string tmpstr;

  if (!m_r_diffs.empty()) {
    take_sample();
  }
  build_index_maps();
  write_header(os, reactions.size());

  for (const auto& sample : m_samples) {
    const auto sim_time = std::get<0>(sample); // time of the reaction
    const s_sample_t& s_sample = std::get<1>(sample);
    const r_sample_t& r_sample = std::get<2>(sample);
    count_species(s_sample, species);
    count_reactions(r_sample, reactions);

    print_stats(sim_time, species, reactions, tmpstr, os);
  }

  return os;
}

void Samples::write(const std::string filename)
{
  std::ofstream ofs;
  ofs.open(filename, std::ofstream::out);
  if (!ofs) return;
  write(ofs);
  ofs.close();
}

/**@}*/
} // end of namespace wcs
