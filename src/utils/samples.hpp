/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef	 __WCS_UTILS_SAMPLE_HPP__
#define	 __WCS_UTILS_SAMPLE_HPP__
#include <string>
#include <unordered_map>
#include <list>
#include <tuple>
#include <fstream>
#include <iostream>
#include "wcs_types.hpp"
#include "reaction_network/network.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

class Samples {

public:
  /// The type of BGL vertex descriptor for graph_t
  using v_desc_t = wcs::Network::v_desc_t;
  using r_prop_t = wcs::Reaction<wcs::Network::v_desc_t>;
  using s_prop_t = wcs::Species;
  using r_desc_t = v_desc_t;
  using s_desc_t = v_desc_t;

  using s_diff_t = wcs::species_cnt_diff_t;
  using r_cnt_t = wcs::species_cnt_t;
  using s_track_t = typename std::pair<s_desc_t, s_diff_t>;
  using r_track_t = typename std::pair<r_desc_t, r_cnt_t>;
  using s_sample_t = std::vector<s_track_t>;
  using r_sample_t = std::vector<r_track_t>;
  using sample_t = std::tuple<sim_time_t, s_sample_t, r_sample_t>;
  using samples_t = std::list<sample_t>;

  Samples();

  void set_time_interval(const sim_time_t t_interval,
                         const sim_time_t t_start = static_cast<sim_time_t>(0));
  void set_iter_interval(const sim_iter_t i_interval,
                         const sim_iter_t i_start = static_cast<sim_iter_t>(0u));

  void record_initial_condition(const std::shared_ptr<wcs::Network>& net_ptr);
  void record_reaction(const sim_time_t t, const r_desc_t r);
  std::ostream& write(std::ostream& os);
  void write(const std::string filename);

protected:
  using s_map_t = typename std::unordered_map<s_desc_t, s_diff_t>;
  using r_map_t = typename std::unordered_map<r_desc_t, r_cnt_t>;

  void take_sample();
  void build_index_maps();
  std::ostream& write_header(std::ostream& os,
                             const size_t num_reactions) const;
  void count_species(const s_sample_t& ss,
                     std::vector<species_cnt_t>& speices) const;
  void count_reactions(const r_sample_t& ss,
                       std::vector<r_cnt_t>& reactions) const;
  std::ostream& print_stats(const sim_time_t sim_time,
                            const std::vector<species_cnt_t>& species,
                            const std::vector<r_cnt_t>& reactions,
                            std::string& tmpstr, std::ostream& os) const;

protected:
  /// Initial species population
  std::vector<species_cnt_t> m_initial_counts;
  /** The pointer to the reaction network being monitored.
   *  Make sure the network object does not get destroyed
   *  while the trace refers to it.
   */
  std::shared_ptr<const wcs::Network> m_net_ptr;
  /// Map a BGL vertex descriptor to the species index
  std::unordered_map<s_desc_t, size_t> m_s_id_map;
  /// Map a BGL vertex descriptor to the reaction index
  std::unordered_map<r_desc_t, size_t> m_r_id_map;
  /**
    * Temporary data structure to keep track of the species count differences
    * over the current sampling interval
    */
  s_map_t m_s_diffs;
  /**
    * Temporary data structure to keep track of the reaction count differences
    * over the current sampling interval
    */
  r_map_t m_r_diffs;

  /// List of samples
  samples_t m_samples;

  sim_iter_t m_start_iter; ///< Simulation iteration of the first record

  sim_iter_t m_cur_iter; ///< Simulation iteration of the current record
  sim_time_t m_cur_time; ///< Simulation time of the current record

  sim_iter_t m_sample_iter_interval; ///< Sampling interval in num of iterations
  sim_time_t m_sample_time_interval; ///< Sampling interval in simulation time

  sim_iter_t m_next_sample_iter; ///< Next iteration to sample
  sim_time_t m_next_sample_time; ///< Next time to sample
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_UTILS_SAMPLE_HPP__
