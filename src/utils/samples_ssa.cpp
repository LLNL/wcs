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

#if defined(WCS_HAS_CEREAL)
#include <cereal/archives/binary.hpp>
#include <cereal/types/list.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/tuple.hpp>
#endif // WCS_HAS_CEREAL

#include <fstream>
#include "utils/samples_ssa.hpp"
#include "utils/to_string.hpp"
#include "utils/exception.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

SamplesSSA::SamplesSSA()
: Trajectory(),
  m_start_iter(static_cast<sim_iter_t>(0u)),
  m_cur_iter(static_cast<sim_iter_t>(0u)),
  m_cur_time(static_cast<sim_time_t>(0)),
  m_sample_iter_interval(static_cast<sim_iter_t>(0u)),
  m_sample_time_interval(static_cast<sim_time_t>(0)),
  m_next_sample_iter(std::numeric_limits<sim_iter_t>::max()),
  m_next_sample_time(std::numeric_limits<sim_time_t>::infinity())
{
}

SamplesSSA::~SamplesSSA()
{}

void SamplesSSA::set_time_interval(const sim_time_t t_interval,
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

void SamplesSSA::set_iter_interval(const sim_iter_t i_interval,
                                const sim_iter_t i_start)
{
  m_sample_iter_interval = i_interval;
  m_start_iter = m_cur_iter = i_start;
  m_next_sample_iter = m_cur_iter + i_interval;
}

/**
 * Build the map from a vertex descriptor to an index of the
 * vector for species and reaction information respectively.
 */
void SamplesSSA::build_index_maps()
{
  Trajectory::build_index_maps();

  if (m_r_id_map.empty()) {
    r_idx_t idx = static_cast<r_idx_t>(0u);
    m_r_id_map.reserve(m_net_ptr->reaction_list().size());

    for (const auto& rd : m_net_ptr->reaction_list()) {
      m_r_id_map[rd] = idx++;
    }
  }
  m_reaction_counts.resize(m_net_ptr->reaction_list().size());
}

/**
 * Accumulate the number of reactions over a certain duration
 * in terms of the simulation time or the number of steps.
 * At the end of duration, sample the state of species counts.
 */

void SamplesSSA::record_step(const sim_time_t t, const SamplesSSA::r_desc_t r)
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

void SamplesSSA::take_sample()
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

 #if defined(WCS_HAS_CEREAL)
  if (++m_cur_record_in_frag >= m_frag_size) {
    flush();
  }
 #endif // WCS_HAS_CEREAL
}

void SamplesSSA::finalize()
{
  if (m_outfile_stem.empty()) {
    m_num_steps = m_samples.size();
    write_header(std::cout);
    write(std::cout);
  } else if ((m_frag_size == static_cast<frag_size_t>(0u)) ||
             (m_frag_size == std::numeric_limits<frag_size_t>::max())) {
    std::ofstream ofs;
    ofs.open((m_outfile_stem + m_outfile_ext), std::ofstream::out);
    if (!ofs) return;

    m_num_steps = m_samples.size();
    write_header(ofs);
    write(ofs);
    ofs.close();
  } else {
  #if defined(WCS_HAS_CEREAL)
    if (!m_r_diffs.empty()) {
      take_sample();
    }

    if (!m_samples.empty()) {
      flush();
    }

    std::ofstream ofs;
    ofs.open((m_outfile_stem + m_outfile_ext), std::ofstream::out);
    write_header(ofs);

    for (frag_id_t i = 0u; i < m_cur_frag_id; ++i) {
      const auto freg_file = m_outfile_stem + '.'
                           + std::to_string(i) + ".cereal";
      std::ifstream is(freg_file, std::ios::binary);
      cereal::BinaryInputArchive archive(is);
      m_samples.clear();
      archive(m_samples);
      write(ofs);
    }

    ofs.close();
  #endif // WCS_HAS_CEREAL
  }
}

size_t SamplesSSA::estimate_tmpstr_size() const
{
  return m_species_counts.size()*cnt_digits +
         m_reaction_counts.size()*cnt_digits;
}

/**
 * Write the header (species labels), and write the initial species population.
 */
std::ostream& SamplesSSA::write_header(std::ostream& os) const
{ // write the header to show the species labels and the initial population
  const auto num_species = m_species_counts.size();
  const auto num_reactions = m_reaction_counts.size();

  sim_iter_t m_num_events = m_cur_iter - m_start_iter;

  std::string ostr = "num_species = " + std::to_string(num_species)
                   + "\tnum_reactions = " + std::to_string(num_reactions)
                   + "\tnum_events = " + std::to_string(m_num_events)
                   + "\nTime: ";

  ostr.reserve(ostr.size() + estimate_tmpstr_size() + 1);

  ostr += m_net_ptr->show_species_labels("")
        + m_net_ptr->show_reaction_labels("")
        + '\n' + std::to_string(0.000000);

  // write the initial population
  for (const auto scnt : m_species_counts) {
    ostr += '\t' + std::to_string(scnt);
  }
  // write the initial reaction distribution
  for (size_t i = 0ul; i < num_reactions; i++) {
    ostr += "\t0";
  }
  os << ostr << '\n';
  return os;
}

void SamplesSSA::count_species(const SamplesSSA::s_sample_t& ss)
{
  for (const auto& s: ss) {
    m_species_counts.at(m_s_id_map.at(std::get<0>(s))) += std::get<1>(s);
  }
}

void SamplesSSA::count_reactions(const SamplesSSA::r_sample_t& rs)
{
  for (const auto& r: rs) {
    m_reaction_counts.at(m_r_id_map.at(std::get<0>(r))) += std::get<1>(r);
  }
}

std::ostream& SamplesSSA::print_stats(
  const sim_time_t sim_time,
  std::string& tmpstr,
  std::ostream& os) const
{
  tmpstr = to_string_in_scientific(sim_time);
  const size_t tstr_sz = tmpstr.size();
  tmpstr.reserve(tstr_sz + 2 + estimate_tmpstr_size());

  for (const auto scnt : m_species_counts) {
    tmpstr += '\t' + std::to_string(scnt);
  }
  for (const auto rcnt : m_reaction_counts) {
    tmpstr += '\t' + std::to_string(rcnt);
  }
  os << tmpstr << '\n';
  return os;
}

std::ostream& SamplesSSA::write(std::ostream& os)
{
  std::string tmpstr;

  if (!m_r_diffs.empty()) {
    take_sample();
  }

  for (const auto& sample : m_samples) {
    const auto sim_time = std::get<0>(sample); // time of the reaction
    const s_sample_t& s_sample = std::get<1>(sample);
    const r_sample_t& r_sample = std::get<2>(sample);
    count_species(s_sample);
    count_reactions(r_sample);

    print_stats(sim_time, tmpstr, os);
  }

  return os;
}

void SamplesSSA::flush()
{
 #if defined(WCS_HAS_CEREAL)
  const auto freg_file = m_outfile_stem + '.' + std::to_string(m_cur_frag_id++)
                       + ".cereal";
  {
    std::ofstream os(freg_file, std::ios::binary);
    cereal::BinaryOutputArchive archive(os);

    archive(m_samples);

    m_cur_record_in_frag = static_cast<frag_size_t>(0u);
  }
  m_num_steps += m_samples.size();
  m_samples.clear();
 #endif // WCS_HAS_CEREAL
}

/**@}*/
} // end of namespace wcs
