#include "reaction_network/network.hpp"
#include "utils/graph_factory.hpp"
#include <type_traits> // is_same<>
#include <algorithm> // lexicographical_compare(), sort()

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  *  @{ */

void Network::load(const std::string graphml_filename)
{
  ::wcs::GraphFactory gfactory;

  if (!gfactory.read_graphml(graphml_filename)) {
    WCS_THROW("Failed to read " + graphml_filename);
    return;
  }
  gfactory.copy_to(m_graph);
}


void Network::init()
{
  const size_t num_vertices = get_num_vertices();

  m_reactions.reserve(num_vertices);
  m_species.reserve(num_vertices);

  v_iter_t vi, vi_end;

  for (boost::tie(vi, vi_end) = boost::vertices(m_graph); vi != vi_end; ++vi) {
    const v_prop_t& v = m_graph[*vi];
    const auto vt = static_cast<v_prop_t::vertex_type>(v.get_typeid());

    if (vt == v_prop_t::_species_) {
      m_species.emplace_back(*vi);
    } else {
      using directed_category = boost::graph_traits<graph_t>::directed_category;
      constexpr bool is_bidirectional 
        = std::is_same<directed_category, boost::bidirectional_tag>::value;

      const v_desc_t reaction = *vi;
      s_involved_t involved_species;

      if constexpr (is_bidirectional) {
        for(const auto ei_in :
            boost::make_iterator_range(boost::in_edges(reaction, m_graph))) {
          v_desc_t reactant = boost::source(ei_in, m_graph);
          involved_species.insert(std::make_pair(m_graph[reactant].get_label(), reactant));
        }
      }

      for(const auto ei_out :
          boost::make_iterator_range(boost::out_edges(reaction, m_graph))) {
        v_desc_t product = boost::target(ei_out, m_graph);
        involved_species.insert(std::make_pair(m_graph[product].get_label(), product));
      }

      m_reactions.emplace_back(reaction);

      auto& r = m_graph[*vi].checked_property< Reaction<v_desc_t> >();
      r.set_rate_inputs(involved_species);
      set_reaction_rate(*vi);
    }
  }

  sort_species();
}

reaction_rate_t Network::set_reaction_rate(const Network::v_desc_t r) const
{
  auto & rprop = m_graph[r].checked_property< Reaction<v_desc_t> >();
  const auto& ri = rprop.get_rate_inputs();
  std::vector<reaction_rate_t> params;
  params.reserve(ri.size()+1u); // reserve space for species count and rate constant
  for (auto vd : ri) { // add species counts here and the rate constant will be appended later
    const auto& s = m_graph[vd].checked_property<Species>();
    params.push_back(static_cast<reaction_rate_t>(s.get_count()));
  }
  return rprop.calc_rate(params);
}

void Network::sort_species()
{
  const auto& g = m_graph;
  std::sort(m_species.begin(), m_species.end(),
            [&g](const Network::v_desc_t lhs, const Network::v_desc_t rhs) {
              const auto lstr = g[lhs].get_label();
              const auto rstr = g[rhs].get_label();
              return std::lexicographical_compare(lstr.begin(), lstr.end(),
                                                  rstr.begin(), rstr.end());
            });
}

Network::v_desc_t Network::find_species(const std::string& label)
{
  const auto& g = m_graph;
  species_list_t::const_iterator it
    = std::lower_bound(m_species.begin(), m_species.end(), label,
            [&g](const Network::v_desc_t lhs, const std::string& l) {
              const auto n = g[lhs].get_label();
              return std::lexicographical_compare(n.begin(), n.end(),
                                                  l.begin(), l.end());
            });
  return *it;
}


size_t Network::get_num_vertices() const
{
  return m_graph.m_vertices.size();
}

size_t Network::get_num_species() const
{
  return m_species.size();
}

size_t Network::get_num_reactions() const
{
  return (get_num_vertices() - get_num_species());
}

size_t Network::get_num_vertices_of_type(const Network::vertex_type vt) const
{
  switch(vt) {
    case v_prop_t::_undefined_: return get_num_vertices();
    case v_prop_t::_species_: return get_num_species();
    case v_prop_t::_reaction_: return get_num_reactions();
  }
  return 0ul;
}

const Network::graph_t& Network::graph() const
{
  return m_graph;
}

const Network::reaction_list_t& Network::reaction_list() const
{
  return m_reactions;
}

const Network::species_list_t& Network::species_list() const
{
  return m_species;
}

/**@}*/
} // end of namespace wcs
