/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __WCS_REACTION_NETWORK_VERTEX_FLAT_HPP__
#define __WCS_REACTION_NETWORK_VERTEX_FLAT_HPP__
// Whole Cell Model Simulator
/** @file
 * \ingroup wcs_reaction_network
 * \brief vertex definition for the reaction graph representation
 */
/** \ingroup wcs_reaction_network
 * \class wcs::VertexFlat
 * \brief vertex definition for the reaction graph representation
 *
 * The bundled vertex property data structure used to define the BGL graph,
 * which represents a reaction network. This in particular provides flat
 * data structure that contains all the member variables that can
 * be used by the derived classes of vertex_base as well as those of
 * vertex_base. This class itself does not have any derived class.
 * It serves as an intermdiate data structure of the bundled property
 * for vertex in parsing GraphML files.
 *
 * \author Jae-Seung Yeom <yeom2@llnl.gov>
 * \date 2019
 */

#include <string>
#include <map>
#include <vector>
#include "wcs_types.hpp"
#include <iostream>

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  @{ */

class VertexFlat {
 public:
  enum vertex_type { _undefined_=0, _species_, _reaction_, _num_vertex_types_ };
  static std::map<vertex_type, std::string> vt_str;
  static std::map<std::string, vertex_type> str_vt;
  /// An alternative to the species_cnt_t that is compatible with graphml
  using s_cnt_t = int;
  /// An alternative type to the reaction_rate_t that is compatible with graphml
  using r_rate_t = double;

  VertexFlat(); ///< BGL requires it to be default constructable

  void set_type(const vertex_type);
  void set_type();
  vertex_type get_type() const;
  int get_typeid() const;
  std::string get_type_str() const;
  void set_label(const std::string& lb);
  std::string get_label() const;
  void set_partition(const partition_id_t pid);
  partition_id_t get_partition() const;

  bool inc_count();
  bool dec_count();
  bool set_count(s_cnt_t c);
  s_cnt_t get_count() const;

  void set_rate_constant(r_rate_t h);
  r_rate_t get_rate_constant() const;
  r_rate_t compute_rate(const std::vector<s_cnt_t>& params);
  r_rate_t get_rate() const;
  void set_rate_formula(const std::string& f);
  const std::string& get_rate_formula() const;

 protected:
  vertex_type m_type;
  int m_typeid;
  std::string m_label; ///< label
  partition_id_t m_pid; ///< The id of the partition to which this edge belongs

  int m_count; ///< copy number of the species
  r_rate_t m_rate; ///< reaction rate
  r_rate_t m_rate_const; ///< rate constant
  std::string m_rate_formula; ///< reaction rate formula

 friend ::wcs::GraphFactory;

 template <typename G>
 friend std::ostream& ::wcs::write_graphviz_of_any_vertex_list(std::ostream& os, const G& g);
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_REACTION_NETWORK_VERTEX_FLAT_HPP__
