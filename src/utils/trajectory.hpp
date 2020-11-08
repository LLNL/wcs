/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef	 __WCS_UTILS_TRAJECTORY_HPP__
#define	 __WCS_UTILS_TRAJECTORY_HPP__
#include <string>
#include <iostream>
#include "sim_methods/update.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

/**
 * Record the sequence of operations performed to show the trajectory of the
 * species population change. Upon completion of simulation, write it into a
 * file. It can shows the population at every change made by each discrete
 * event in the order (tracing), or only show the aggregate change at a regular
 * interval (sampling).
 * During simulation, trajectory data are kept in a memory buffer. The buffer
 * can be configured to flush out to temprary files, called fragments, as it
 * fills up to an amount predefined by users.
 * At finalization, the history of population change is reconstructed from the
 * buffer or the fragment files, and written into a final trajectory file.
 */
class Trajectory {
public:
  /// The type of BGL vertex descriptor for graph_t
  using s_prop_t = wcs::Species;
  using r_desc_t = wcs::Network::v_desc_t;
  using r_prop_t = wcs::Network::r_prop_t;
  using r_cnt_t = wcs::species_cnt_t;
  using frag_id_t = size_t;
  using map_desc2idx_t = wcs::Network::map_desc2idx_t;

  Trajectory(const std::shared_ptr<wcs::Network>& net_ptr);
  Trajectory(const Trajectory& other) = default;
  Trajectory(Trajectory&& other) = default;
  Trajectory& operator=(const Trajectory& other) = default;
  Trajectory& operator=(Trajectory&& other) = default;

  virtual ~Trajectory();
  void set_outfile(const std::string outfile = "",
                   const frag_size_t frag_size = default_frag_size);

  virtual void initialize();
  virtual void record_step(const sim_time_t t, const r_desc_t r);
  virtual void record_step(const sim_time_t t, cnt_updates_t&& updates);
  virtual void record_step(const sim_time_t t, conc_updates_t&& updates);
  virtual void finalize(const sim_time_t t) = 0;

protected:
  void record_initial_condition();
  virtual std::ostream& write_header(std::ostream& os) const = 0;
  virtual std::ostream& write(std::ostream& os) = 0;
  virtual void flush();

protected:
  /// Initial species population
  std::vector<species_cnt_t> m_species_counts;
  /** The pointer to the reaction network being monitored.
   *  Make sure the network object does not get destroyed
   *  while the trace refers to it.
   */
  std::shared_ptr<const wcs::Network> m_net_ptr;
  /// Map a BGL vertex descriptor to the reaction index
  const map_desc2idx_t* m_r_id_map;
  /// Map a BGL vertex descriptor to the species index
  const map_desc2idx_t* m_s_id_map;

  /// Output file name stem (the part without extention)
  std::string m_outfile_stem;
  /// Output file name extension
  std::string m_outfile_ext;

  /// Number of records per output fragment
  frag_size_t m_frag_size;
  /// The id of the current fragment
  frag_id_t m_cur_frag_id;
  /// The id of the current record in current fragment
  frag_size_t m_cur_record_in_frag;
  /// Total number of steps recorded
  size_t m_num_steps;
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_UTILS_TRAJECTORY_HPP__
