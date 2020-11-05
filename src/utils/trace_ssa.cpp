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
#endif // WCS_HAS_CEREAL

#include <fstream>
#include "utils/trace_ssa.hpp"
#include "utils/exception.hpp"
#include "utils/to_string.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

TraceSSA::TraceSSA(const std::shared_ptr<wcs::Network>& net_ptr)
: Trajectory(net_ptr)
{}

TraceSSA::~TraceSSA()
{}

void TraceSSA::record_step(const sim_time_t t, const r_desc_t r)
{
  m_trace.emplace_back(std::make_pair(t, m_net_ptr->reaction_d2i(r)));

 #if defined(WCS_HAS_CEREAL)
  if (++m_cur_record_in_frag >= m_frag_size) {
    flush();
  }
 #endif // WCS_HAS_CEREAL
}

void TraceSSA::initialize()
{
  Trajectory::initialize();
  m_reaction_counts.resize(m_net_ptr->get_num_reactions());
}

void TraceSSA::finalize()
{
  if (m_outfile_stem.empty()) {
    m_num_steps = m_trace.size();
    write_header(std::cout);
    write(std::cout);
  } else if ((m_frag_size == static_cast<frag_size_t>(0u)) ||
             (m_frag_size == std::numeric_limits<frag_size_t>::max())) {
    std::ofstream ofs;
    ofs.open((m_outfile_stem + m_outfile_ext), std::ofstream::out);
    if (!ofs) return;

    m_num_steps = m_trace.size();
    write_header(ofs);
    write(ofs);
    ofs.close();
  } else {
  #if defined(WCS_HAS_CEREAL)
    if (!m_trace.empty()) {
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
      m_trace.clear();
      archive(m_trace);
      write(ofs);
    }

    ofs.close();
  #endif // WCS_HAS_CEREAL
  }
}

size_t TraceSSA::estimate_tmpstr_size() const
{
  return m_species_counts.size()*cnt_digits +
         m_reaction_counts.size()*cnt_digits;
}

/**
 * Write the header (species labels), and write the initial species population.
 */
std::ostream& TraceSSA::write_header(std::ostream& os) const
{ // write the header to show the species labels and the initial population
  const auto num_species = m_net_ptr->get_num_species();
  const auto num_reactions = m_net_ptr->get_num_reactions();

  std::string ostr = "num_species = " + std::to_string(num_species)
                   + "\tnum_reactions = " + std::to_string(num_reactions)
                   + "\tnum_events = " + std::to_string(m_num_steps)
                   + "\nTime: ";

  ostr += m_net_ptr->show_species_labels("")
        + "\tReaction"
        + m_net_ptr->show_reaction_labels("")
        + '\n' + std::to_string(0.000000);

  ostr.reserve(ostr.size() + estimate_tmpstr_size() + 4);

  // write the initial population
  for (const auto scnt : m_species_counts) {
    ostr += '\t' + std::to_string(scnt);
  }
  ostr += "\tNA";
  // write the initial reaction distribution
  for (const auto rcnt : m_reaction_counts) {
    ostr += '\t' + std::to_string(rcnt);
  }
  os << ostr << '\n';
  return os;
}

std::ostream& TraceSSA::print_stats(const sim_time_t sim_time,
                                    const std::string rlabel,
                                    std::string& tmpstr,
                                    std::ostream& os) const
{
  tmpstr = to_string_in_scientific(sim_time);
  const size_t tstr_sz = tmpstr.size();
  tmpstr.reserve(tstr_sz + 2 + estimate_tmpstr_size());

  for (const auto scnt : m_species_counts) {
    tmpstr += '\t' + std::to_string(scnt);
  }
  tmpstr += '\t' + rlabel;
  for (const auto rcnt : m_reaction_counts) {
    tmpstr += '\t' + std::to_string(rcnt);
  }
  os << tmpstr << '\n';
  return os;
}

std::ostream& TraceSSA::write(std::ostream& os)
{
  if (!m_net_ptr) {
    WCS_THROW("Invaid pointer for reaction network.");
  }
  const wcs::Network::graph_t& g = m_net_ptr->graph();

  std::string tmpstr;

  // species population state
  auto& species = m_species_counts;

  trace_t::const_iterator it = m_trace.begin();
  trace_t::const_iterator it_end = m_trace.cend();

  const auto& r_desc_map = m_net_ptr->reaction_list();

  for (; it != it_end; it++) {
    const auto sim_time = it->first; // time of the reaction
    // BGL vertex descriptor of the reaction
    r_desc_t vd_reaction = r_desc_map.at(it->second);
    m_reaction_counts.at(m_r_id_map->at(vd_reaction)) ++;

    // product species
    for (const auto ei_out : boost::make_iterator_range(boost::out_edges(vd_reaction, g))) {
      const auto vd_product = boost::target(ei_out, g);
      if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
        // in case that there are other type of vertices than species or reaction
        if (g[vd_product].get_type() != wcs::Vertex::_species_) continue;
      }
      const auto stoichio = g[ei_out].get_stoichiometry_ratio();
      species.at(m_s_id_map->at(vd_product)) += stoichio;
    }

    // reactant species
    for (const auto ei_in : boost::make_iterator_range(boost::in_edges(vd_reaction, g))) {
      const auto vd_reactant = boost::source(ei_in, g);
      if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
        // in case that there are other type of vertices than species or reaction
        if (g[vd_reactant].get_type() != wcs::Vertex::_species_) continue;
      }
      const auto stoichio = g[ei_in].get_stoichiometry_ratio();
      species.at(m_s_id_map->at(vd_reactant)) -= stoichio;
    }
    print_stats(sim_time, g[vd_reaction].get_label(), tmpstr, os);
  }

  return os;
}

void TraceSSA::flush()
{
 #if defined(WCS_HAS_CEREAL)
  const auto freg_file = m_outfile_stem + '.' + std::to_string(m_cur_frag_id++)
                       + ".cereal";
  {
    std::ofstream os(freg_file, std::ios::binary);
    cereal::BinaryOutputArchive archive(os);

    archive(m_trace);

    m_cur_record_in_frag = static_cast<frag_size_t>(0u);
  }
  m_num_steps += m_trace.size();
  m_trace.clear();
 #endif // WCS_HAS_CEREAL
}

/**@}*/
} // end of namespace wcs
