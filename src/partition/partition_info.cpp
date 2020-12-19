/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include <type_traits>

#include <partition/partition_info.hpp>
#include <utils/detect_methods.hpp>
#include <utils/exception.hpp>
#include <utils/write_graphviz.hpp>

namespace wcs {
/** \addtogroup wcs_partition
 *  *  @{ */

Partition_Info::Partition_Info(const std::shared_ptr<const wcs::Network>& net,
                               const size_t num_partitions,
                               const load_t min_load)
: m_net(net), m_num_parts(num_partitions), m_min_load(min_load)
{
  if constexpr (!is_bidirectional) {
    WCS_THROW("Bi-directional BGL network is required.");
    return;
  }

  if (m_num_parts == 0ul) {
    reset_num_partitions();
  }

  m_num_local.resize(n_types, p_sum_t(m_num_parts));

  m_sum_indegree.resize(n_types, p_sum_t(m_num_parts));
  m_sum_outdegree.resize(n_types, p_sum_t(m_num_parts));

  m_max_indegree.resize(n_types, p_max_t(m_num_parts));
  m_max_outdegree.resize(n_types, p_max_t(m_num_parts));
  m_max_degree.resize(n_types, p_max_t(m_num_parts));

  m_load.resize(n_types, p_load_t(m_num_parts, load_t{0}));
  m_comm_out.resize(n_types, part_adj_t(m_num_parts));
  m_comm_in.resize(n_types, part_adj_t(m_num_parts));
}

void Partition_Info::reset_num_partitions()
{
  if (!m_net) {
    return;
  }

  const graph_t& g = m_net->graph();
  v_iter_t vi, vi_end;
  boost::tie(vi, vi_end) = boost::vertices(g);
  auto& num_parts = m_num_parts;

  std::function<void (const v_desc_t&)> max_part_id
    = [&] (const v_desc_t& v) {
    const auto pid = static_cast<size_t>(g[v].get_partition());

    if (num_parts <= pid) {
      num_parts = pid + 1;
    }
  };

  std::for_each(vi, vi_end, max_part_id);
}

void Partition_Info::scan(bool verbose)
{
  const graph_t& g = m_net->graph();

  if (m_net->get_partition_id() == unassigned_partition) {
    using std::operator<<;
    std::cerr << "This network is not partitioned." << std::endl;
    return;
  }

  { // local vertices, degrees, and load
    size_t v_idx = 0ul;
    std::function<void (const v_desc_t&)> part_info_vertex
      = [&, this] (const v_desc_t& v)
    {
      const auto& vp = g[v];
      // vertex type
      const auto vt = vp.get_type();
      // map the vertex type to a vector index
      const auto type = static_cast<size_t>(vt-1);
      // partition id
      const auto pid = vp.get_partition();

      const auto num_in_edges = boost::in_degree(v, g);
      const auto num_out_edges = boost::out_degree(v, g);
      const auto degree = num_in_edges + num_out_edges;

      if (verbose) {
        using std::operator<<;
        std::cout << " V " << v_idx++ << "\t"
                  << vp.get_label() << "\t"
                  << "\tP" << pid
                  << '\t' << num_in_edges
                  << ' ' << num_out_edges
                  << std::endl;
      }
      m_num_local.at(type).at(pid) ++;
      m_sum_indegree.at(type).at(pid) += num_in_edges;
      m_sum_outdegree.at(type).at(pid) += num_out_edges;

      auto& max_in_deg = m_max_indegree.at(type).at(pid);
      auto& max_out_deg = m_max_outdegree.at(type).at(pid);
      auto& max_deg = m_max_degree.at(type).at(pid);

      if (max_in_deg.first == num_in_edges) {
        max_in_deg.second.insert(v);
      } else if (max_in_deg.first < num_in_edges) {
        max_in_deg = std::make_pair(num_in_edges, std::set{v});
      }

      if (max_out_deg.first == num_out_edges) {
        max_out_deg.second.insert(v);
      } else if (max_out_deg.first < num_out_edges) {
        max_out_deg = std::make_pair(num_out_edges, std::set{v});
      }

      if (max_deg.first == degree) {
        max_deg.second.insert(v);
      } else if (max_deg.first < degree) {
        max_deg = std::make_pair(degree, std::set{v});
      }

      if (vt == v_prop_t::_reaction_) {
        const auto& rp = vp.property<r_prop_t>();
        const auto load = rp.get_rate();
        m_load.at(type).at(pid) += ((load == 0)? m_min_load : load);
      } else {
        m_load.at(type).at(pid) += m_min_load;
      }
    };

    v_iter_t vi, vi_end;
    boost::tie(vi, vi_end) = boost::vertices(g);

    std::for_each(vi, vi_end, part_info_vertex);
  }

  { // edges and communication volume
    std::function<void (const e_desc_t&)> edge_writer
      = [&, this] (const e_desc_t& e)
    {
      const auto src = boost::source(e, g);
      const auto& vp_src = g[src];
      const auto pid_src = vp_src.get_partition();

      const auto dst = boost::target(e, g);
      const auto& vp_dst = g[dst];
      const auto pid_dst = vp_dst.get_partition();

      if (pid_src == pid_dst) {
        return;
      }

      const auto vt_src = vp_src.get_type();
      const auto vt_dst = vp_dst.get_type();

      if constexpr (wcs::Vertex::_num_vertex_types_  > 3) {
        // Assuming the reaction vertex is the only one with compute load
        if (vt_src != v_prop_t::_reaction_ && vt_dst != v_prop_t::_reaction_) {
          return;
        }
      }

      const auto& rp_src = vp_src.property<r_prop_t>();
      const auto& rp_dst = vp_dst.property<r_prop_t>();
      auto rate = reaction_rate_t{0};

      if (vt_src == v_prop_t::_reaction_) {
        const auto deg = boost::out_degree(src, g);
        rate = rp_src.get_rate() / deg;
      } else {
        const auto deg = boost::in_degree(dst, g);
        rate = rp_dst.get_rate() / deg;
      }
      const auto type_src = static_cast<size_t>(vt_src-1);
      const auto type_dst = static_cast<size_t>(vt_dst-1);
      m_comm_out.at(type_src).at(pid_src)[pid_dst] += rate;
      m_comm_in.at(type_dst).at(pid_dst)[pid_src] += rate;
    };

    e_iter_t ei, ei_end;
    boost::tie(ei, ei_end) = boost::edges(g);

    std::for_each(ei, ei_end, edge_writer);
  }

  return;
}

template<size_t NumVTypes, typename SumPerPart>
static void show_sum_per_part(
  const SumPerPart& spp,
  const std::string& title,
  size_t nparts)
{
  using std::operator<<;
  using namespace std;
  using value_t = typename SumPerPart::value_type::value_type;

  cout << title << endl;
  for (size_t p = 0ul; p < nparts; ++p) {
    auto total = value_t{0};
    string str = 'P' + to_string(p);
    for (size_t t = 0ul; t < NumVTypes; ++t) {
      const auto n = spp.at(t).at(p);
      total += n;
      str += "\t" + Vertex::vt_str.at(static_cast<Vertex::vertex_type>(t+1))
                  + ' ' + to_string(n);
    }
    cout << (str + "\tTot " + to_string(total)) << endl;
  }
}

template<size_t NumVTypes, typename MaxPerPart>
static void show_max_per_part(
  const MaxPerPart& mpp,
  const std::string& title,
  size_t nparts,
  const Partition_Info::graph_t& g,
  bool brief_max)
{
  using std::operator<<;
  using namespace std;

  cout << title << endl;
  for (size_t p = 0ul; p < nparts; ++p) {
    string str = 'P' + to_string(p) + '\t';
    for (size_t t = 0ul; t < NumVTypes; ++t) {
      const auto mx = mpp.at(t).at(p);
      if (mx.first <= 0ul) {
        continue;
      }
      const auto vt_str
        = Vertex::vt_str.at(static_cast<Vertex::vertex_type>(t+1));
      if (brief_max) {
        str += " <" + vt_str + ' ' + to_string(mx.first) + '>';
      } else {
        for (const auto& vd: mx.second) {
          const auto v_label = g[vd].get_label();
          str += ' ' + vt_str
               + " <" + v_label + ", " + to_string(mx.first) + ">";
        }
      }
    }
    cout << str << endl;
  }
}

template<size_t NumVTypes, typename CommVolPerPart>
static void show_comm_volume_per_part(
  const CommVolPerPart& cvpp,
  const std::string& title,
  const bool direction, // outgoing = true, incoming = false
  size_t nparts)
{
  using std::operator<<;
  using namespace std;
  using value_t = typename CommVolPerPart::value_type::value_type::mapped_type;

  cout << title << endl;
  for (size_t p_src = 0ul; p_src < nparts; ++p_src) {
    std::vector<value_t> total (nparts, static_cast<value_t>(0));
    string str = "  P" + to_string(p_src) + (direction? " ->" : " <-");
    for (size_t p_dst = 0ul; p_dst < nparts; ++p_dst) {
      for (size_t t = 0ul; t < NumVTypes; ++t) {
        const auto& adj = cvpp.at(t).at(p_src);
        Partition_Info::part_edges_t::const_iterator it = adj.find(p_dst);
        if (it != adj.cend()) {
          total[p_dst] += it->second;
        }
      }
      if (total[p_dst] > value_t{0}) {
        str += " <P" + to_string(p_dst) + ", " + to_string(total[p_dst]) + ">";
      }
    }
    cout << str << endl;
  }
}

void Partition_Info::report(bool brief_max) const
{
  using std::operator<<;
  using namespace std;
  const graph_t& g = m_net->graph();

  cout << std::endl << "Num parts: " << m_num_parts << endl;

  std::cout << "======================================";
  show_sum_per_part<n_types>(m_num_local,
                             "\nNum of local vertices:",
                             m_num_parts);

  show_sum_per_part<n_types>(m_sum_indegree,
                             "\nSum of indegrees of local vertices:",
                             m_num_parts);

  show_sum_per_part<n_types>(m_sum_outdegree,
                             "\nSum of outdegrees of local vertices:",
                             m_num_parts);

  std::cout << "--------------------------------------";
  show_max_per_part<n_types>(m_max_indegree,
                             "\nMaximum of indegrees of local vertices:",
                             m_num_parts, g, brief_max);

  show_max_per_part<n_types>(m_max_outdegree,
                             "\nMaximum of outdegrees of local vertices:",
                             m_num_parts, g, brief_max);

  show_max_per_part<n_types>(m_max_degree,
                             "\nMaximum of degrees of local vertices:",
                             m_num_parts, g, brief_max);

  std::cout << "--------------------------------------";
  show_sum_per_part<n_types>(m_load,
                             "\nCompute load:",
                             m_num_parts);

  show_comm_volume_per_part<n_types>(m_comm_out,
                             "\nCommunication volume (outgoing):",
                             true,
                             m_num_parts);

  show_comm_volume_per_part<n_types>(m_comm_in,
                             "\nCommunication volume (incoming):",
                             false,
                             m_num_parts);
  std::cout << "======================================" << std::endl << std::endl;
}

/**@}*/
} // end of namespace wcs
