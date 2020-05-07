/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __WCS_REACTION_NETWORK_VERTEX_HPP__
#define __WCS_REACTION_NETWORK_VERTEX_HPP__
// Whole Cell Model Simulator
/** @file
 * \ingroup wcs_reaction_network
 * \brief vertex definition for the reaction graph representation
 */
/** \ingroup wcs_reaction_network
 * \class wcs::Vertex
 * \brief vertex definition for the reaction graph representation
 *
 * The bundled vertex property data structure used to define the BGL graph,
 * which represents a reaction network.
 *
 * \author Jae-Seung Yeom <yeom2@llnl.gov>
 * \date 2019
 */

#include <string>
#include <map>
#include <iostream>
#include "wcs_types.hpp"
#include "bgl.hpp"
#include "reaction_network/vertex_flat.hpp"
#include "reaction_network/vertex_property_base.hpp"
#include "reaction_network/species.hpp"
#include "reaction_network/reaction.hpp"

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  *  @{ */

class Vertex {
 public:
  enum vertex_type { _undefined_=0, _species_, _reaction_, _num_vertex_types_ };
  static std::map<vertex_type, std::string> vt_str;
  static std::map<std::string, vertex_type> str_vt;

  Vertex(); ///< BGL requires it to be default constructable
  Vertex(const Vertex& rhs);
  Vertex(Vertex&& rhs) noexcept;
  Vertex& operator=(const Vertex& rhs);
  Vertex& operator=(Vertex&& rhs) noexcept;
  template <typename G> Vertex(const VertexFlat& rhs, const G& g);
  virtual ~Vertex();
  std::unique_ptr<Vertex> clone() const;

  vertex_type get_type() const;
  int get_typeid() const;
  std::string get_type_str() const;
  void set_label(const std::string& lb);
  std::string get_label() const;

  template <typename P> P& property() const;
  template <typename P> P& checked_property() const;

 protected:
  void set_type(const vertex_type);
  void set_type();
  void reset(Vertex& obj);

 private:
  virtual Vertex* clone_impl() const;

 protected:
  vertex_type m_type; ///< The vertex type
  int m_typeid; ///< The vertex type in integer form used for GraphML parsing
  std::string m_label; ///< The vertex label

  /**
   * The pointer to the detailed property object, which is polymorphic.
   * This allows to use different property bundle structure types for reaction
   * and species vertices even when the BGL does not support polymorphic vertex
   * property bundle type.
   */
  std::unique_ptr<VertexPropertyBase> m_p;

 friend ::wcs::GraphFactory;

 template <typename G>
 friend std::ostream& ::wcs::write_graphviz_of_any_vertex_list(std::ostream& os, const G& g);
};


template <typename G>
Vertex::Vertex(const VertexFlat& flat, const G& g)
: m_type(static_cast<vertex_type>(flat.get_typeid())),
  m_typeid(flat.get_typeid()),
  m_label(flat.get_label()),
  m_p(nullptr)
{
  m_typeid = static_cast<int>(m_type);

  switch(m_type) {
    case _species_: {
        m_p = std::unique_ptr<Species>(new Species);
        dynamic_cast<Species*>(m_p.get())->set_count(flat.get_count());
        break;
      }
    case _reaction_: {
        using v_desc_t = typename boost::graph_traits<G>::vertex_descriptor;
        m_p = std::unique_ptr< Reaction<v_desc_t> >(new Reaction<v_desc_t>);
        dynamic_cast<Reaction<v_desc_t>*>(m_p.get())->set_rate_constant(flat.get_rate_constant());
        dynamic_cast<Reaction<v_desc_t>*>(m_p.get())->set_rate_formula(flat.get_rate_formula());
        break;
      }
    default:
      break;
  }
}

template <typename P> P& Vertex::property() const
{
  return *dynamic_cast<P*>(m_p.get());
}

template <typename P> P& Vertex::checked_property() const
{
  auto ptr = dynamic_cast<P*>(m_p.get());
  if (ptr == nullptr) {
    WCS_THROW("Attempted to dereference a wrong type of property porinter.");
  }
  return *ptr;
}

/**@}*/
} // end of namespace wcs
#endif // __WCS_REACTION_NETWORK_VERTEX_HPP__
