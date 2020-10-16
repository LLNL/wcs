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
#include "utils/trace_generic.hpp"
#include "utils/exception.hpp"
#include "utils/to_string.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

TraceGeneric::~TraceGeneric()
{}

void TraceGeneric::record_step(const sim_time_t t, cnt_updates_t&& updates)
{
  m_trace.emplace_back(std::make_pair(t, std::forward<cnt_updates_t>(updates)));

 #if defined(WCS_HAS_CEREAL)
  if (++m_cur_record_in_frag >= m_frag_size) {
    flush();
  }
 #endif // WCS_HAS_CEREAL
}

/**
 * Build the map from a vertex descriptor to an index of the
 * vector for species respectively.
 */
void TraceGeneric::build_index_maps()
{
  if (m_s_id_map.empty()) {
    r_idx_t idx = static_cast<r_idx_t>(0u);
    m_s_id_map.reserve(m_net_ptr->species_list().size());

    for (const auto& sd : m_net_ptr->species_list()) {
      m_s_id_map[sd] = idx++;
    }
  }
}

void TraceGeneric::finalize()
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

/**
 * Write the header (species labels), and write the initial species population.
 */
std::ostream& TraceGeneric::write_header(std::ostream& os) const
{ // write the header to show the species labels and the initial population
  const auto num_species = m_net_ptr->get_num_species();
  const auto num_reactions = m_net_ptr->get_num_reactions();

  std::string ostr = "num_species = " + std::to_string(num_species)
                   + "\tnum_reactions = " + std::to_string(num_reactions)
                   + "\tnum_events = " + std::to_string(m_num_steps)
                   + "\nTime: ";

  ostr.reserve(ostr.size() + m_net_ptr->get_num_species()*cnt_digits + 4);

  ostr += m_net_ptr->show_species_labels("")
        + "\tOperation\n" + std::to_string(0.000000);

  // write the initial population
  for (const auto scnt : m_species_counts) {
    ostr += '\t' + std::to_string(scnt);
  }
  os << ostr << "\tNA\n";
  return os;
}

size_t TraceGeneric::estimate_tmpstr_size() const
{
  return m_species_counts.size()*cnt_digits;
}

std::ostream& TraceGeneric::print_stats(const sim_time_t sim_time,
                                    const std::string elabel,
                                    std::string& tmpstr,
                                    std::ostream& os) const
{
  tmpstr = to_string_in_scientific(sim_time);
  const size_t tstr_sz = tmpstr.size();
  tmpstr.reserve(tstr_sz + 2 + estimate_tmpstr_size());

  for (const auto scnt : m_species_counts) {
    tmpstr += '\t' + std::to_string(scnt);
  }
  os << tmpstr << '\t' + elabel + '\n';

  return os;
}

std::ostream& TraceGeneric::write(std::ostream& os)
{
  if (!m_net_ptr) {
    WCS_THROW("Invaid pointer for reaction network.");
  }

  std::string tmpstr;

  // species population state
  auto& species = m_species_counts;

  trace_t::const_iterator it = m_trace.begin();
  trace_t::const_iterator it_end = m_trace.cend();

  for (; it != it_end; it++) {
    const auto sim_time = it->first; // time of event
    const auto& updates = it->second; // population updates

    for (const auto& u : updates) {
      species.at(m_s_id_map.at(u.first))
        += static_cast<species_cnt_diff_t>(u.second);
    }

    print_stats(sim_time, "", tmpstr, os);
  }

  return os;
}

void TraceGeneric::flush()
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
