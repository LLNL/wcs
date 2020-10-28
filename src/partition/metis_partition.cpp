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

// Only enable when METIS is available
#if defined(WCS_HAS_METIS)
#include <set>
#include <cstddef> // NULL used in Metis
#include <string>
#include "partition/metis_partition.hpp"
#include "utils/exception.hpp"

namespace wcs {
/** \addtogroup wcs_partition
 *  @{ */

Metis_Partition::Metis_Partition(const Metis_Params& mp)
: m_p(mp),
  m_num_edges(0ul),
  m_num_vertices(0ul) {}

size_t Metis_Partition::get_num_edges() const
{
  return m_num_edges;
}

size_t Metis_Partition::get_num_vertices() const
{
  return m_num_vertices;
}


void Metis_Partition::prepare()
{
  if (!m_p.m_rnet) {
    WCS_THROW("Invalid pointer to reaction network!");
    return;
  }
  m_num_vertices = (m_p.m_rnet)->get_num_vertices();
  build_map_from_desc_to_idx();
  populate_adjacny_list();
  populate_vertex_info();
}


void Metis_Partition::build_map_from_desc_to_idx()
{
  const graph_t& graph = (m_p.m_rnet)->graph();

  m_vd2idx.clear();
  m_vd2idx.reserve(get_num_vertices());
  m_idx2vd.clear();

  m_num_edges = 0ul;
  idx_t idx = static_cast<idx_t>(0);

  v_iter_t vi, vi_end;

  for (boost::tie(vi, vi_end) = boost::vertices(graph); vi != vi_end; ++vi) {
    const v_desc_t& vertex = *vi;
    m_idx2vd.emplace(idx, vertex);
    m_vd2idx.emplace(vertex, idx++);
    //if constexpr (is_bidirectional) {
    //  m_num_edges += boost::in_degree(vertex, graph);
    //}
    m_num_edges += boost::out_degree(vertex, graph);
  }
  //if constexpr (is_bidirectional) {
  //  m_num_edges =  m_num_edges / 2;
  //}
}


void Metis_Partition::populate_adjacny_list()
{
  const graph_t& graph = (m_p.m_rnet)->graph();

  // Offset array of adjacency list in CSR
  m_xadj.clear();
  m_xadj.reserve(get_num_vertices()+1);
  m_xadj.push_back(static_cast<idx_t>(0));

  // Ajacency list in CSR
  m_adjncy.clear();
  m_adjncy.reserve(get_num_edges()*2);

  v_iter_t vi, vi_end;

  // Populate m_xadj and m_adjncy
  if constexpr (is_bidirectional) {
    for (boost::tie(vi, vi_end) = boost::vertices(graph); vi != vi_end; ++vi) {
      std::vector<idx_t> neighbors;

      const v_desc_t& vertex = *vi;
      //const wcs::Network::v_prop_t& v = graph[vertex];
      for(const auto ei_in :
          boost::make_iterator_range(boost::in_edges(vertex, graph)))
      {
        neighbors.push_back(m_vd2idx.at(boost::source(ei_in, graph)));
      }

      for(const auto ei_out :
          boost::make_iterator_range(boost::out_edges(vertex, graph)))
      {
        neighbors.push_back(m_vd2idx.at(boost::target(ei_out, graph)));
      }
      std::sort(neighbors.begin(), neighbors.end());
      m_xadj.push_back(m_xadj.back() + neighbors.size());
      m_adjncy.insert(m_adjncy.end(), neighbors.begin(), neighbors.end());
    }
  } else {
    std::vector< std::set<idx_t> > neighbors_list(get_num_vertices());

    for (boost::tie(vi, vi_end) = boost::vertices(graph); vi != vi_end; ++vi) {
      const v_desc_t& vertex = *vi;
      //const wcs::Network::v_prop_t& v = graph[vertex];

      for(const auto ei_out :
          boost::make_iterator_range(boost::out_edges(vertex, graph)))
      {
        const auto v = m_vd2idx.at(vertex);
        const auto n = m_vd2idx.at(boost::target(ei_out, graph));
        std::set<idx_t>& v_neighbors = neighbors_list.at(v);
        std::set<idx_t>& n_neighbors = neighbors_list.at(n);
        v_neighbors.insert(n);
        n_neighbors.insert(v);
      }
    }

    for (size_t i = 0ul; i < neighbors_list.size(); ++i) {
      const auto& neighbors = neighbors_list[i];
      m_xadj.push_back(m_xadj.back() + neighbors.size());
      m_adjncy.insert(m_adjncy.end(), neighbors.begin(), neighbors.end());
    }
  }
}


/*
 * Vertex sizes are used in computing communication volume in Metis. This method
 * uses the same value as the vertex weight for it. The rationale behind it is
 * that the cost of each reaction event is roughly the same, and the
 * computational cost of a reaction vertex is proportional to the the number of
 * the event or the rate of the reaction of the vertex. Similarly, the cost of
 * communication for each event message is the same. Thus, the total cost of
 * communication is proportional to the number of messages or the rate of the
 * reaction.
 *
 * This method only considers single-valueid vertex weight (load balancing
 * constraint). In case of a species vertex, we assign the minimum value
 * possible possible for its weight, which is 1 as it needs to be a non-zero
 * positive integral value. In case of a reaction vertex, we linearly convert
 * the floating-point rate into an integer value.
 */
void Metis_Partition::populate_vertex_info()
{
  const graph_t& graph = (m_p.m_rnet)->graph();
  const auto n_vertices = get_num_vertices();
  // Minimum vertext weight (needs to be non-zero positive interger)
  const idx_t vwgt_min = m_p.m_vwgt_min;
  const idx_t vwgt_max = m_p.m_vwgt_max;

  if (n_vertices  == 0ul) {
    WCS_THROW("Empty graph!");
    return;
  }

  // Vertex weights
  m_vwgt.clear();
  m_vwgt.reserve(n_vertices);

  // Vertex sizes (used in com_puting communication volumne)
  m_vsize.clear();

  if (vwgt_max == vwgt_min) {
    return;
  }

  // Min, max, and sum of reaction rates
  const auto [r_min, r_max, r_sum] = m_p.m_rnet->find_min_max_rate();

  // Upper bound of vertex weight representation
  constexpr idx_t vwgt_ub = std::numeric_limits<idx_t>::max() - 1;

  if (r_min < static_cast<reaction_rate_t>(0)) {
    WCS_THROW("Minimum reaction rate must be a non-negative number!");
    return;
  }
  if (r_max <= static_cast<reaction_rate_t>(0)) {
    WCS_THROW("Maximum reaction rate must be a positive number!");
    return;
  }
  if (r_sum <= static_cast<reaction_rate_t>(0)) {
    WCS_THROW("Sum of reaction rate must be a positive number!");
    return;
  }

  const reaction_rate_t r_ratio = r_max / r_sum;

  if ((vwgt_ub < n_vertices) ||
      (static_cast<double>(vwgt_ub)/n_vertices < static_cast<double>(vwgt_min))) {
    WCS_THROW("(Num_vertices X min_vertex_weight) is too large!");
    return;
  }

  v_iter_t vi, vi_end;

  /* In case of a reaction vertex, we convert the floating-point rate `r'
   * into an integer value `vwgt' given `width' and `vwgt_min' as
   * `vwgt = vwgt_min + width * (r/r_max)'
   *
   * `width' is the difference between the largest vertext weight and the
   * minimum (vwgt_min). It should satisfy (sum(vwgt) < vwgt_ub) where
   * `sum(vwgt) = vwgt_min * n_vertices + width * (r_sum/r_max)'
   *
   * Given vwgt_min, the largest width can be computed as follows.
   * Suppose `sum(vwgt)' is the maximum allowed value `vwgt_ub'. Then,
   * `width = (vwgt_ub - vwgt_min * n_vertices) * (r_max / r_sum)'
   *
   * Non-reaction vertex will be assigned `vwgt_min'.
   */

  // Largest possible width
  idx_t width = (vwgt_ub - vwgt_min * n_vertices) * r_ratio;
  // Adjust width if requested
  if ((vwgt_max > static_cast<idx_t>(0)) &&
      (width + vwgt_min > vwgt_max))
  {
    width = vwgt_max - vwgt_min;
  }

  //const idx_t vwgt_sum
  //  = vwgt_min * n_vertices + static_cast<idx_t>(width * (r_sum/r_max));

  using std::operator<<;
  for (boost::tie(vi, vi_end) = boost::vertices(graph); vi != vi_end; ++vi) {
    const v_prop_t& v = graph[*vi];
    const auto vt = static_cast<v_prop_t::vertex_type>(v.get_typeid());

    if (vt == v_prop_t::_reaction_) {
      const auto& rv = graph[*vi]; // vertex (property) of the reaction
      const auto& rp = rv.property<r_prop_t>(); // detailed vertex property data
      const auto r = rp.get_rate();
      const idx_t vwgt = vwgt_min + static_cast<idx_t>(width * r/r_max);
      m_vwgt.push_back(vwgt);
    } else {
      m_vwgt.push_back(vwgt_min);
    }
  }

  if (m_p.m_ratio_w2s > 0.0) {
    m_vsize = m_vwgt;

    if (m_p.m_ratio_w2s != 1.0) {
      const double ratio_w2s = m_p.m_ratio_w2s;
      std::for_each(m_vsize.begin(), m_vsize.end(),
        [ratio_w2s](idx_t& s) {
          s = (s <= ratio_w2s)?
                static_cast<idx_t>(1) :
                static_cast<idx_t>(s/ratio_w2s);
        });
    }
  }
}


bool Metis_Partition::check_run(const int ret, const bool verbose)
{
  std::string msg;

  switch (ret) {
    case METIS_OK:
      msg = "METIS_OK";
      break;
    case METIS_ERROR_INPUT:
      msg = "METIS_ERROR_INPUT";
      break;
    case METIS_ERROR_MEMORY:
      msg = "METIS_ERROR_MEMORY";
      break;
    case METIS_ERROR:
      msg = "METIS_ERROR";
      break;
  }

  if (verbose) {
  using std::operator<<;
    std::cout << msg << std::endl;
  }
  return (ret == METIS_OK);
}


bool Metis_Partition::run(std::vector<idx_t>& parts, idx_t& objval)
{ // Unfortunately, the Metis function used below does not take const parameter
  idx_t nvtxs = static_cast<idx_t>(get_num_vertices());

  parts.clear();
  parts.resize(nvtxs, static_cast<idx_t>(0));
  objval = static_cast<idx_t>(0);

  idx_t* vwgt_ptr = (m_vwgt.empty()? NULL : m_vwgt.data());
  idx_t* vsize_ptr = (m_vsize.empty()? NULL : m_vsize.data());

  int ret
    = METIS_PartGraphKway(&nvtxs, &(m_p.m_nvwghts), m_xadj.data(), m_adjncy.data(),
                          vwgt_ptr, vsize_ptr, NULL, &(m_p.m_nparts), NULL,
                          NULL, m_p.m_opts.data(), &objval, parts.data());

  return check_run(ret);
}


void Metis_Partition::print_params() const
{
  using std::operator<<;
  std::cout << "Num vertices: " << get_num_vertices() << std::endl;
  std::cout << "Num edges: " << get_num_edges() << std::endl;
  m_p.print();
}


void Metis_Partition::print_metis_inputs(
  const std::vector<idx_t>& xadj, const std::vector<idx_t>& adjncy,
  const std::vector<idx_t>& vwgt, const std::vector<idx_t>& vsize,
  std::ostream& os)
{
  using std::operator<<;
  os << std::endl << "xadj:" << std::endl << "\t";
  for (const auto x : xadj) {
    os << ' ' << x;
  }

  os << std::endl << "adjncy:" << std::endl << "\t";
  for (const auto a : adjncy) {
    os << ' ' << a;
  }

  os << std::endl << "vwgt:" << std::endl << "\t";
  for (const auto w : vwgt) {
    os << ' ' << w;
  }

  os << std::endl << "vsize:" << std::endl << "\t";
  for (const auto s : vsize) {
    os << ' ' << s;
  }
  os << std::endl;
}


void Metis_Partition::print_metis_graph(std::ostream& os) const {
  print_metis_inputs(m_xadj, m_adjncy,
                     m_vwgt, m_vsize, os);
}


void Metis_Partition::print_adjacency(std::ostream& os) const
{
  using std::operator<<;
  const graph_t& graph = (m_p.m_rnet)->graph();

  os << std::endl << "Adjacency:";
  for (size_t i = 1; i < m_xadj.size(); ++i) {
    os << std::endl;
    const auto& vd = m_idx2vd.at(i-1);
    os << '\t' << graph[vd].get_label() << " ->";
    for (auto j = m_xadj[i-1]; j < m_xadj[i]; ++j) {
      const auto& vd_connected = m_idx2vd.at(m_adjncy[j]);
      os << ' ' << graph[vd_connected].get_label();
    }
  }

  os << std::endl;
}

const std::unordered_map<Network::v_desc_t, idx_t>&
Metis_Partition::get_map_from_desc_to_idx() const
{
  return m_vd2idx;
}

const std::map<idx_t, Network::v_desc_t>&
Metis_Partition::get_map_from_idx_to_desc() const
{
  return m_idx2vd;
}

/**@}*/
} // end of namespace wcs
#endif // defined(WCS_HAS_CONFIG)
