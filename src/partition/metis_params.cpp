#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

// Only enable when METIS is available
#if defined(WCS_HAS_METIS)
#include <iostream>
#include <map>
#include <string>
#include "partition/metis_params.hpp"
#include "utils/exception.hpp"

namespace wcs {
/** \addtogroup wcs_partition
 *  @{ */

Metis_Params::Metis_Params()
: m_nvwghts(0), m_nparts(0),
  m_vwgt_min(static_cast<idx_t>(1)),
  m_vwgt_max(m_vwgt_min),
  m_ratio_w2s(0.0),
  m_rnet(nullptr)
{
  METIS_SetDefaultOptions(m_opts.data());

  m_opts[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY; // Partitioning Schemes
  //m_opts[METIS_OPTION_CTYPE] = METIS_CTYPE_RM; // Coarsening by random matching
  m_opts[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM; // Coarsening by sorted heavy-edge matching

  m_opts[METIS_OPTION_NITER] = 10; // Number of refinement iterations
  m_opts[METIS_OPTION_NCUTS] = 1; // Number of different partitionings to compute
  m_opts[METIS_OPTION_SEED] = 7177;
  m_opts[METIS_OPTION_MINCONN] = 0; // Minimize the maximum degree of subdomain graph
  m_opts[METIS_OPTION_UFACTOR] = 300; // Maximum allowed load imbalance
  m_opts[METIS_OPTION_NUMBERING] = 0; // C-style numbering that starts from 0

  make_options_consistent();
}

void Metis_Params::make_options_consistent()
{
  if (m_ratio_w2s == 0.0) {
    // Vertex size is not used. Thus, the partitioning objective is edge-cut
    m_opts[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
  } else {
    // Partitioning objective: comm volumne
    m_opts[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
  }
  if (m_vwgt_max == m_vwgt_min) {
    m_nvwghts = 0;
  } else {
    m_nvwghts = 1;
  }
}

void Metis_Params::limit_max_vertex_weight(idx_t w)
{
  if (w < static_cast<idx_t>(0)) {
    WCS_THROW("Upper-bound of vertex weight must be non-negative!");
    return;
  } else {
    m_vwgt_max = w;
  }
  make_options_consistent();
}

void Metis_Params::set_ratio_of_vertex_weight_to_size(double r)
{
  if (r < 0.0) {
    WCS_THROW("Ratio of vertex weight to size must be non-negative!");
    return;
  }

  m_ratio_w2s = r;
  make_options_consistent();
}

/// Specify the number of desired partitions and the input graph
bool Metis_Params::set(idx_t np, std::shared_ptr<wcs::Network> rnet)
{
  m_rnet = rnet;
  if (!m_rnet) {
    WCS_THROW("Invalid reaction network pointer!");
    return false;
  }
  m_nparts = np;
  return true;
}

/// Set relevant Metis partitioning options
bool Metis_Params::set_options(mobjtype_et objective, mctype_et coarsening,
                               idx_t niter, idx_t seed, bool minconn, idx_t ufactor,
                               uint8_t dbglvl)
{
  using std::operator<<;
#if 0 // Enable if the parameters are not the enum types defined in Metis
  bool ok = (objective == METIS_OBJTYPE_CUT || // by edge cut
             objective == METIS_OBJTYPE_VOL || // by comm volume
             objective == METIS_OBJTYPE_NODE); // by node
  if (!ok) {
    WCS_THROW("Invalid option for partitioning objective!");
    return false;
  }

  // Random matching or sorted heavy-edge matching
  ok = (coarsening == METIS_CTYPE_RM || coarsening == METIS_CTYPE_SHEM);
  if (!ok) {
    WCS_THROW("Invalid option for coarsening!");
    return false;
  }
#endif

  if (niter < static_cast<idx_t>(1)) {
    WCS_THROW("Invalid number of refinement iterations (" \
              + std::to_string(niter) + ")!");
    return false;
  }

  if ((ufactor <= static_cast<idx_t>(0)) &&
      (ufactor > static_cast<idx_t>(1000))) {
    WCS_THROW("Invalid option for ufactor (" \
              + std::to_string(ufactor) + ")!");
    return false;
  }

  // Partitioning objective
  m_opts[METIS_OPTION_OBJTYPE] = static_cast<idx_t>(objective);
  // Coarsening
  m_opts[METIS_OPTION_CTYPE] = static_cast<idx_t>(coarsening);
  // Number of refinement iterations
  m_opts[METIS_OPTION_NITER] = static_cast<idx_t>(niter);
  // Set seed for the randomization RNG
  m_opts[METIS_OPTION_SEED] = seed;
  // Minimize the maximum degree of subdomain graph
  m_opts[METIS_OPTION_MINCONN] = static_cast<idx_t>(minconn);
  // Maximum allowed load imbalance
  m_opts[METIS_OPTION_UFACTOR] = static_cast<idx_t>(ufactor);
  m_opts[METIS_OPTION_DBGLVL] = static_cast<idx_t>(dbglvl);

  make_options_consistent();

  return true;
}

idx_t Metis_Params::get_seed() const
{
  return (m_opts[METIS_OPTION_SEED]);
}

void Metis_Params::print() const
{
  static const std::map <idx_t, const std::string> objective_str = {
    {static_cast<idx_t>(METIS_OBJTYPE_CUT), "METIS_OBJTYPE_CUT"},
    {static_cast<idx_t>(METIS_OBJTYPE_VOL), "METIS_OBJTYPE_VOL"},
    {static_cast<idx_t>(METIS_OBJTYPE_NODE), "METIS_OBJTYPE_NODE"}
  };

  static const std::map <idx_t, const std::string> coarsening_str =  {
    {static_cast<idx_t>(METIS_CTYPE_RM), "METIS_CTYPE_RM"},
    {static_cast<idx_t>(METIS_CTYPE_SHEM), "METIS_CTYPE_SHEM"}
  };

  using std::operator<<;

  std::cout << "------ Metis params set ------" << std::endl;
  std::cout << " - Num vertex weights (LB constraints): " << m_nvwghts << std::endl;
  std::cout << " - Num partitions: " << m_nparts << std::endl;
  std::cout << " - Vertex weight lower-bound: " << m_vwgt_min << std::endl;
  std::cout << " - Vertex weight upper-bound: " << m_vwgt_max << std::endl;
  std::cout << " - Ratio of vertex weight to size: " << m_ratio_w2s << std::endl;

  std::cout << " - Objective: " << objective_str.at(m_opts[METIS_OPTION_OBJTYPE]) << std::endl;
  std::cout << " - Coarsening: " << coarsening_str.at(m_opts[METIS_OPTION_CTYPE]) << std::endl;
  std::cout << " - Num refinement iterations: " << m_opts[METIS_OPTION_NITER] << std::endl;
  std::cout << " - Metis RN seed: " << get_seed() << std::endl;
  std::cout << " - Minconn: " << m_opts[METIS_OPTION_MINCONN] << std::endl;
  std::cout << " - Ufactor: " << m_opts[METIS_OPTION_UFACTOR] << std::endl;
}

/**@}*/
} // end of namespace wcs
#endif // defined(WCS_HAS_METIS)
