#include <iostream>
#include <boost/filesystem.hpp>
#include <utils/graph_factory.hpp>

extern "C" {
#if HAVE_CONFIG_H
#include "config.h"
#endif
}

namespace wcs {
/** \addtogroup wcs_utils
 *  *  @{ */

void GraphFactory::setup_dynamic_property ()
{
  m_dp.property("v_label", get(&v_prop_t::m_label, m_g));
  m_dp.property("v_type",  get(&v_prop_t::m_typeid, m_g));
  m_dp.property("s_count", get(&v_prop_t::m_count, m_g));
  m_dp.property("r_const", get(&v_prop_t::m_rate_const, m_g));
  m_dp.property("r_rate",  get(&v_prop_t::m_rate_formula, m_g));
  m_dp.property("e_label", get(&e_prop_t::m_label, m_g));
  m_dp.property("e_stoic", get(&e_prop_t::m_stoichio, m_g));
}


GraphFactory::GraphFactory ()
{
  setup_dynamic_property ();
}

GraphFactory::GraphFactory (const GraphFactory &o)
{
  m_g = o.m_g;
  m_dp = o.m_dp;
}

const GraphFactory::graph_t &GraphFactory::graph () const
{
  return m_g;
};


/*! Load reaction network from a graphml file
 *
 *  \param ifn   filename of the input reaction network in graphml
 *  \return      true on success; false otherwise
 */
bool GraphFactory::read_graphml(const std::string &ifn)
{
  bool ok = true;
  std::ifstream in_file(ifn.c_str());
  if (!in_file.good ())
    return false;

  try {
    boost::read_graphml(in_file, m_g, m_dp);
  } catch (boost::graph_exception &e) {
    std::cerr << e.what () << std::endl;
    ok = false;
  }

  in_file.close ();

  // Set m_type of each vertex according to m_typeid.
  // This is because graphml does not support enum type
  // and thus we rely on the integer representation of
  // the m_type, which is m_typeid.
  std::function<void (const v_desc_t&)> set_vertex_type
    = [&] (const v_desc_t& v) {
      m_g[v].set_type();
    };

  typename boost::graph_traits<graph_t>::vertex_iterator vi, vi_end;
  boost::tie(vi, vi_end) = boost::vertices(m_g);
  std::for_each(vi, vi_end, set_vertex_type);

  return ok;
}

/**@}*/
} // end of namespace wcs
