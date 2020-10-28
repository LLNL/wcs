/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef  __WCS_METIS_PARAMS__
#define  __WCS_METIS_PARAMS__

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#if defined(WCS_HAS_METIS)
#include <metis.h>
#include <array>
#include <memory>
#include "reaction_network/network.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

struct Metis_Params {
  idx_t m_nvwghts; ///< Number of vertex weights (constraints)
  idx_t m_nparts; ///< Number of partitions
  /// Lower-bound on the vertext weight (needs to be non-zero positive interger)
  idx_t m_vwgt_min;
  idx_t m_vwgt_max; ///< Upper-bound on the vertex weight
  double m_ratio_w2s; ///< ratio of vertex weight to vertex size
  std::shared_ptr<const wcs::Network> m_rnet; ///< Reaction network

  std::array<idx_t, METIS_NOPTIONS> m_opts; ///< Metis options

  Metis_Params();
  Metis_Params(const Metis_Params& other) = default;
  Metis_Params(Metis_Params&& other) = default;
  Metis_Params& operator=(const Metis_Params& other) = default;
  Metis_Params& operator=(Metis_Params&& other) = default;

  /// Set essential parameters. i.e., # partitions, and the graph to partition
  bool set(idx_t np, std::shared_ptr<const wcs::Network> rnet);
  /// Set options for Metis partitioner
  bool set_options(mobjtype_et objective, mctype_et coarsening, idx_t niter,
                   idx_t seed, bool minconn, idx_t ufactor, uint8_t dbglvl);
  void make_options_consistent();

  /**
   * Specify how large an integer vertex weight can be when mapped from the
   * floating-point reaction rate. When zero is given as the parameter, the
   * largest possible value of the upper-bound of the vertex weight is used.
   * The upper-bound will be used to populate the vertex weight list (the
   * load balance constraint input). A reaction rate maps to a number between
   * the lower-bound and the upper-bound.
   *
   * The default of the upper-bound is the same value as the lower-bound.
   * In such a case, the vertex weight is not used in partitioning.
   */
  void limit_max_vertex_weight(idx_t w);
  /**
   * The vertex size is used to compute the volumne of the communication out of
   * the vertex. In our reaction network problem, it is set as the weight of
   * the vertex divided by the ratio given here. If zero is given (which is
   * the default), then the vertex size is not used in partitioning.
   */ 
  void set_ratio_of_vertex_weight_to_size(double r);
  
  void print() const;
};

/**@}*/
} // end of namespace wcs
#endif // defined(WCS_HAS_METIS)
#endif //  __WCS_METIS_PARAMS__
