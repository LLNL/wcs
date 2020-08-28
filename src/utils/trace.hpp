/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef	 __WCS_UTILS_TRACE_HPP__
#define	 __WCS_UTILS_TRACE_HPP__
#include <map>
#include <list>
#include <string>
#include <fstream>
#include <iostream>
#include "wcs_types.hpp"
#include "reaction_network/network.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

class Trace {
public:
  /// The type of BGL vertex descriptor for graph_t
  using v_desc_t = wcs::Network::v_desc_t;
  using r_prop_t = wcs::Reaction<wcs::Network::v_desc_t>;
  using s_prop_t = wcs::Species;
  using r_desc_t = v_desc_t;
  using s_desc_t = v_desc_t;
  using tentry_t = typename std::pair<sim_time_t, r_desc_t>;
  using trace_t = std::list<tentry_t>;

  virtual ~Trace();
  void record_reaction(const sim_time_t t, const r_desc_t r);
  void record_initial_condition(const std::shared_ptr<wcs::Network>& net_ptr);
  void pop_back();

  virtual void write(const std::string filename);
  virtual std::ostream& write(std::ostream& os);

protected:
  virtual void build_index_maps();
  virtual std::ostream& write_header(std::ostream& os) const;
  virtual void count_reaction(r_desc_t r);
  virtual size_t estimate_tmpstr_size() const;
  virtual std::ostream& print_stats(const sim_time_t sim_time,
                                    const std::vector<species_cnt_t>& species,
                                    const std::string rlabel,
                                    std::string& tmpstr, std::ostream& os) const;

protected:
  /// Trace records
  trace_t m_trace;
  /// Initial species population
  std::vector<species_cnt_t> m_initial_counts;
  /** The pointer to the reaction network being monitored.
   *  Make sure the network object does not get destroyed
   *  while the trace refers to it.
   */
  std::shared_ptr<const wcs::Network> m_net_ptr;
  /// Map a BGL vertex descriptor to the species index
  std::unordered_map<s_desc_t, size_t> m_s_id_map;
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_UTILS_TRACE_HPP__
