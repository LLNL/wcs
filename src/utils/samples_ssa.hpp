/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef	 __WCS_UTILS_SAMPLE_SSA_HPP__
#define	 __WCS_UTILS_SAMPLE_SSA_HPP__
#include <list>
#include <tuple>
#include <iostream>
#include "utils/trajectory.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

/**
 * Record aggregate changes in the species population at a regular interval
 * predefined by users. The interval can be specified as either the number
 * of operations or the amount of simulation time passed. Upon completion of
 * simulation, write a file that shows how the species population has
 * changed over time at a regular interval.
 * During simulation, operation records are kept in a memory buffer. The
 * buffer can be configured to flush out to a file as it fills up to an
 * amount predefined by users.
 * At finalization, the history of species population change is reconstructed
 * from the buffer or the history fragment files, and written into a file.
 */
class SamplesSSA : public Trajectory {

public:
  using s_diff_t = wcs::species_cnt_diff_t;
  // TODO: convert s_track_t and r_track_t to use integer instead of descriptors
  using s_track_t = typename std::pair<s_desc_t, s_diff_t>;
  using r_track_t = typename std::pair<r_desc_t, r_cnt_t>;
  using s_sample_t = std::vector<s_track_t>;
  using r_sample_t = std::vector<r_track_t>;
  using sample_t = std::tuple<sim_time_t, s_sample_t, r_sample_t>;
  using samples_t = std::list<sample_t>;

  SamplesSSA();
  ~SamplesSSA() override;

  void set_time_interval(const sim_time_t t_interval,
                         const sim_time_t t_start = static_cast<sim_time_t>(0));
  void set_iter_interval(const sim_iter_t i_interval,
                         const sim_iter_t i_start = static_cast<sim_iter_t>(0u));

  void record_step(const sim_time_t t, const r_desc_t r) override;
  void finalize() override;

protected:
  using s_map_t = typename std::unordered_map<s_desc_t, s_diff_t>;
  using r_map_t = typename std::unordered_map<r_desc_t, r_cnt_t>;

  void build_index_maps() override;
  void take_sample();
  size_t estimate_tmpstr_size() const;
  std::ostream& write_header(std::ostream& os) const override;
  void count_species(const s_sample_t& ss);
  void count_reactions(const r_sample_t& ss);
  std::ostream& print_stats(const sim_time_t sim_time,
                            std::string& tmpstr, std::ostream& os) const;
  std::ostream& write(std::ostream& os) override;
  void flush() override;

protected:
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

  /// Show how many times each reaction fires
  std::vector<r_cnt_t> m_reaction_counts;

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
#endif // __WCS_UTILS_SAMPLE_SSA_HPP__
