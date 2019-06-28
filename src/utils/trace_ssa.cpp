#include "utils/trace_ssa.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  *  @{ */

/**
 * Build the map from a vertex descriptor to an index of the
 * vector for species and reaction information respectively.
 */
void TraceSSA::build_index_maps()
{
  Trace::build_index_maps();

  if (m_r_id_map.empty()) {
    size_t idx = 0ul;
    m_r_id_map.reserve(m_net_ptr->reaction_list().size());

    for (const auto& rd : m_net_ptr->reaction_list()) {
      m_r_id_map[rd] = idx++;
    }
  }
  m_reaction_counts.resize(m_net_ptr->reaction_list().size());
}

/**
 * Write the header (species labels), and write the initial species population.
 */
std::ostream& TraceSSA::write_header(std::ostream& os) const
{ // write the header to show the species labels and the initial population
  std::string ostr("Time: ");
  ostr.reserve(ostr.size() + m_net_ptr->get_num_species()*20);

  ostr += m_net_ptr->show_species_labels("")
        + "\tReaction"
        + m_net_ptr->show_reaction_labels("")
        + '\n' + std::to_string(0.000000);

  // write the initial population
  for (const auto scnt : m_initial_counts) {
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

void TraceSSA::count_reaction(r_desc_t r)
{
  m_reaction_counts.at(m_r_id_map.at(r)) ++;
}

size_t TraceSSA::estimate_tmpstr_size() const
{
  return m_initial_counts.size()*10;
}

std::ostream& TraceSSA::print_stats(const sim_time_t sim_time,
                                    const std::vector<species_cnt_t>& species,
                                    const std::string rlabel,
                                    std::string& tmpstr,
                                    std::ostream& os) const
{
  tmpstr = std::to_string(sim_time);
  tmpstr.reserve(estimate_tmpstr_size());

  for (const auto scnt : species) {
    tmpstr += '\t' + std::to_string(scnt);
  }
  tmpstr += '\t' + rlabel;
  for (const auto rcnt : m_reaction_counts) {
    tmpstr += '\t' + std::to_string(rcnt);
  }
  os << tmpstr << '\n';
  return os;
}

/**@}*/
} // end of namespace wcs
