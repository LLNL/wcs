/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __WCS_TYPES_HPP__
#define __WCS_TYPES_HPP__
#include <iostream>
#include <utility> // std::move
#include <memory>  // std::unique_ptr
#include <limits>  // std::numeric_limits

namespace wcs {
/** \addtogroup wcs_global
 *  @{ */

using species_cnt_t = unsigned int;
using species_cnt_diff_t = int;
using reaction_rate_t = double;
using sim_time_t = double;
using sim_iter_t = unsigned;
using stoic_t = int;
using partition_id_t = unsigned;
using concentration_t = double;
using v_idx_t = unsigned; ///< vertex index type

constexpr partition_id_t unassigned_partition
  = std::numeric_limits<partition_id_t>::max();

template <typename G>
std::ostream& write_graphviz_of_any_vertex_list(std::ostream& os, const G& g);

class GraphFactory;

constexpr size_t num_in_edges_to_reserve = 8ul;
/// ceil(log10(2^sizeof(species_cnt_t))) + sizeof('\t')
constexpr size_t cnt_digits = 21ul;

/// The default number of records per tracing/sampling output file fragment
using frag_size_t = unsigned;
constexpr frag_size_t default_frag_size = 16384u;

/**@}*/
} // end of namespace wcs
#endif // __WCS_TYPES_HPP__
