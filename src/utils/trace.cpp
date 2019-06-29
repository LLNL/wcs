#include "utils/trace.hpp"
#include "utils/exception.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  *  @{ */

void Trace::record_reaction(const sim_time_t t, const Trace::r_desc_t r)
{
  m_trace.emplace_back(std::make_pair(t, r));
}

void Trace::record_initial_condition(const std::shared_ptr<wcs::Network>& net_ptr)
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
    m_initial_counts[i++] = sp.get_count();
  }
}

void Trace::pop_back()
{
  m_trace.pop_back();
}

void Trace::write(const std::string filename)
{
  std::ofstream ofs;
  ofs.open(filename, std::ofstream::out);
  if (!ofs) return;
  write(ofs);
  ofs.close();
}

/**
 * Build the map from a vertex descriptor to an index of the
 * vector for species respectively.
 */
void Trace::build_index_maps()
{
  if (m_s_id_map.empty()) {
    size_t idx = 0ul;
    m_s_id_map.reserve(m_net_ptr->species_list().size());

    for (const auto& sd : m_net_ptr->species_list()) {
      m_s_id_map[sd] = idx++;
    }
  }
}

/**
 * Write the header (species labels), and write the initial species population.
 */
std::ostream& Trace::write_header(std::ostream& os) const
{ // write the header to show the species labels and the initial population
  std::string ostr("Time: ");
  ostr.reserve(ostr.size() + m_net_ptr->get_num_species()*20);

  ostr += m_net_ptr->show_species_labels("")
        + "\tReaction\n" + std::to_string(0.000000);

  // write the initial population
  for (const auto scnt : m_initial_counts) {
    ostr += '\t' + std::to_string(scnt);
  }
  os << ostr << "\tNA\n";
  return os;
}

void Trace::count_reaction(r_desc_t r) {}

size_t Trace::estimate_tmpstr_size() const
{
  return m_initial_counts.size()*10;
}

std::ostream& Trace::print_stats(const sim_time_t sim_time,
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
  os << tmpstr << '\t' + rlabel + '\n';

  return os;
}

std::ostream& Trace::write(std::ostream& os)
{
  if (!m_net_ptr) {
    WCS_THROW("Invaid pointer for reaction network.");
  }
  const wcs::Network::graph_t& g = m_net_ptr->graph();

  build_index_maps();
  write_header(os);

  std::string tmpstr;

  // species population state
  std::vector<species_cnt_t> species(m_initial_counts);

  // TODO: capability to specify the window of time preiod
  trace_t::const_iterator it = m_trace.begin();
  trace_t::const_iterator it_end = m_trace.cend();

  sim_time_t sim_time = 0.0; // simulation time

  for (; it != it_end; it++) {
    sim_time += it->first; // time of the reaction
    r_desc_t vd_reaction = it->second; // BGL vertex descriptor of the reaction
    count_reaction(vd_reaction);

    // product species
    for (const auto ei_out : boost::make_iterator_range(boost::out_edges(vd_reaction, g))) {
      const auto vd_product = boost::target(ei_out, g);
      if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
        // in case that there are other type of vertices than species or reaction
        if (g[vd_product].get_type() != wcs::Vertex::_species_) continue;
      }
      const auto stoichio = g[ei_out].get_stoichiometry_ratio();
      species.at(m_s_id_map.at(vd_product)) += stoichio;
    }

    // reactant species
    for (const auto ei_in : boost::make_iterator_range(boost::in_edges(vd_reaction, g))) {
      const auto vd_reactant = boost::source(ei_in, g);
      if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
        // in case that there are other type of vertices than species or reaction
        if (g[vd_reactant].get_type() != wcs::Vertex::_species_) continue;
      }
      const auto stoichio = g[ei_in].get_stoichiometry_ratio();
      species.at(m_s_id_map.at(vd_reactant)) -= stoichio;
    }
    print_stats(sim_time, species, g[vd_reaction].get_label(), tmpstr, os);
  }

  return os;
}

/**@}*/
} // end of namespace wcs
