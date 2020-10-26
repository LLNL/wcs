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

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#include <string>
#include <map>
#include <iostream>
#include "wcs_types.hpp"
#include "bgl.hpp"
#include "reaction_network/vertex_flat.hpp"
#include "reaction_network/vertex_property_base.hpp"
#include "reaction_network/species.hpp"
#include "reaction_network/reaction.hpp"
#include "utils/sbml_utils.hpp"


#if defined(WCS_HAS_SBML)
#include <sbml/SBMLTypes.h>
#include <sbml/common/extern.h>
#endif // defined(WCS_HAS_SBML)

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  @{ */

typedef reaction_rate_t ( * rate_function_pointer)(const std::vector<reaction_rate_t>&);

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

  #if defined(WCS_HAS_SBML)
  template <typename G>
  Vertex(const LIBSBML_CPP_NAMESPACE::Species& species, const G& g);

  template <typename G>
  Vertex(const LIBSBML_CPP_NAMESPACE::Model& model, const
   LIBSBML_CPP_NAMESPACE::Reaction& reaction, const G& g, const
  std::function<reaction_rate_t (const std::vector<reaction_rate_t>&)>&  reaction_function_rate);
  #endif // defined(WCS_HAS_SBML)

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
 friend std::ostream&
 ::wcs::write_graphviz_of_any_vertex_list(std::ostream& os, const G& g);
};


template <typename G>
Vertex::Vertex(const VertexFlat& flat, const G& g)
: m_type(static_cast<vertex_type>(flat.get_typeid())),
  m_typeid(flat.get_typeid()),
  m_label(flat.get_label()),
  m_p(nullptr)
{
  switch(m_type) {
    case _species_: {
        m_p = std::unique_ptr<Species>(new Species);
        dynamic_cast<Species*>(m_p.get())->set_count(flat.get_count());
        break;
      }
    case _reaction_: {
        using v_desc_t = typename boost::graph_traits<G>::vertex_descriptor;
        m_p = std::unique_ptr< Reaction<v_desc_t> >(new Reaction<v_desc_t>);
        dynamic_cast<Reaction<v_desc_t>*>(m_p.get())->
          set_rate_constant(flat.get_rate_constant());
        dynamic_cast<Reaction<v_desc_t>*>(m_p.get())->
          set_rate_formula(flat.get_rate_formula());
    break;
      }
    default:
      break;
  }
}

#if defined(WCS_HAS_SBML)
template <typename G>
Vertex::Vertex(const LIBSBML_CPP_NAMESPACE::Species& species, const G& g)
: m_type(_species_),
  m_typeid(static_cast<int>(_species_)),
  m_label(species.getIdAttribute()),
  m_p(nullptr)
{
  m_p = std::unique_ptr<Species>(new Species);

  if  (!isnan(species.getInitialAmount())) {
    dynamic_cast<Species*>(m_p.get())->
      set_count(static_cast<species_cnt_t>(species.getInitialAmount()));
    // TODO: species.getInitialConcentration() should be multiplied by
    // compartment volume and molarity (depending on the unit of
    // concentration) should be converted to count
  } else  if (!isnan(species.getInitialConcentration())) {
    dynamic_cast<Species*>(m_p.get())->
      set_count(static_cast<species_cnt_t>(species.getInitialConcentration()));
  }
}

template <typename G>
Vertex::Vertex(
  const LIBSBML_CPP_NAMESPACE::Model& model,
  const LIBSBML_CPP_NAMESPACE::Reaction& reaction,
  const G& g,
  const std::function<reaction_rate_t (const std::vector<reaction_rate_t>&)>&  reaction_function_rate
  )
: m_type(_reaction_),
  m_typeid(static_cast<int>(_reaction_)),
  m_label(reaction.getIdAttribute()),
  m_p(nullptr)
{
  using v_desc_t = typename boost::graph_traits<G>::vertex_descriptor;
  m_p = std::unique_ptr< Reaction<v_desc_t> >(new Reaction<v_desc_t>);
  dynamic_cast<Reaction<v_desc_t>*>(m_p.get())->set_rate_constant(1);

  using reaction_parameters = std::unordered_set<std::string>;
  typename reaction_parameters::const_iterator mpit;
  reaction_parameters mpset;
  reaction_parameters pset;

  std::string formula = SBML_formulaToString(reaction.getKineticLaw()->getMath());

  //remove spaces from formula
  formula.erase(remove(formula.begin(), formula.end(), ' '), formula.end());

  std::string wholeformula("");
  //Add parameters
  const LIBSBML_CPP_NAMESPACE::ListOfParameters* parameter_list
    = model.getListOfParameters();
  unsigned int parametersSize = parameter_list->size();

  sbml_utils sbml_o;
  pset=sbml_o.get_reaction_parameters(model, reaction);

  // Create an unordered_set for all model parameters
  for (unsigned int pi = 0u; pi < parametersSize; pi++) {
    const LIBSBML_CPP_NAMESPACE::Parameter* parameter = parameter_list->get(pi);
    mpset.insert(parameter->getIdAttribute());
  }

  //Add parameters
  for  (const std::string& x: pset) {
    std::string s_label = x;
    mpit = mpset.find(s_label) ;
    if (mpit != mpset.end()) {
      std::stringstream ss;
      ss << parameter_list->get(s_label)->getValue();
      std::string parametervalue = ss.str();
      wholeformula = wholeformula + "var " + s_label
                   + " := " + parametervalue +  "; ";
    }
  }

  //Add compartnents
  const LIBSBML_CPP_NAMESPACE::ListOfCompartments* compartment_list
    = model.getListOfCompartments();
  unsigned int compartmentsSize = compartment_list->size();

  for (unsigned int ic = 0u; ic < compartmentsSize; ic++) {
    const LIBSBML_CPP_NAMESPACE::Compartment* compartment
      = compartment_list->get(ic);
    std::string toFindPar(compartment->getIdAttribute());
    size_t posPar = formula.find(toFindPar);

    std::stringstream ss;
    ss << compartment->getSize();
    std::string parametervalue = ss.str();

    if (posPar != std::string::npos) {
      wholeformula = wholeformula + "var " + compartment->getIdAttribute()
                   + " := " + parametervalue +  "; ";
    }
  }

  wholeformula = wholeformula + "m_rate := " + formula + ";";

  dynamic_cast<Reaction<v_desc_t>*>(m_p.get())->set_rate_formula(wholeformula);

  #if !defined(WCS_HAS_EXPRTK)
  dynamic_cast<Reaction<v_desc_t>*>(m_p.get())->ReactionBase::set_calc_rate_fn(reaction_function_rate);
  #endif // !defined(WCS_HAS_EXPRTK)
}
#endif // defined(WCS_HAS_SBML)

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
