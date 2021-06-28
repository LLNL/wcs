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
  Vertex(const LIBSBML_CPP_NAMESPACE::Model& model, const LIBSBML_CPP_NAMESPACE::Species& species, const G& g);

  template <typename G>
  Vertex(const LIBSBML_CPP_NAMESPACE::Model& model, const
    LIBSBML_CPP_NAMESPACE::Reaction& reaction, const G& g, const
    std::function<reaction_rate_t (const std::vector<reaction_rate_t>&)>&
    reaction_function_rate);
  #endif // defined(WCS_HAS_SBML)

  virtual ~Vertex();
  std::unique_ptr<Vertex> clone() const;

  vertex_type get_type() const;
  int get_typeid() const;
  std::string get_type_str() const;
  void set_label(const std::string& lb);
  std::string get_label() const;
  void set_partition(const partition_id_t pid);
  partition_id_t get_partition() const;
  void set_annotation(const std::string& f);
  std::string get_annotation() const;

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
  partition_id_t m_pid; ///< The id of the partition to which this edge belongs
  std::string m_annotation; ///< vertex annotation

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
  m_pid(unassigned_partition),
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
Vertex::Vertex(const LIBSBML_CPP_NAMESPACE::Model& model, 
const LIBSBML_CPP_NAMESPACE::Species& species, const G& g)
: m_type(_species_),
  m_typeid(static_cast<int>(_species_)),
  m_label(species.getIdAttribute()),
  m_pid(unassigned_partition),
  m_p(nullptr)
{
  m_p = std::unique_ptr<Species>(new Species);

  if  (species.isSetInitialAmount()) {
    dynamic_cast<Species*>(m_p.get())->
      set_count(static_cast<species_cnt_t>(species.getInitialAmount()));
    // TODO: if units of species are mole then species Initial Amount has 
    // to be multiplied by Avogadro number

    // TODO: species.getInitialConcentration() should be multiplied by
    // compartment volume and molarity (depending on the unit of
    // concentration) should be converted to count
  } else  if (species.isSetInitialConcentration()) {
    const LIBSBML_CPP_NAMESPACE::ListOfCompartments* compartment_list
    = model.getListOfCompartments();
    const LIBSBML_CPP_NAMESPACE::ListOfUnitDefinitions* unit_definition_list
    = model.getListOfUnitDefinitions(); 
    double compartment_size = compartment_list->get(species.getCompartment())->getSize();
    std::string compartment_unit ;
    if (compartment_list->get(species.getCompartment())->isSetUnits()){
      compartment_unit = compartment_list->get(species.getCompartment())->getUnits();
    } else if (compartment_list->get(species.getCompartment())->getSpatialDimensions() == 3) {
      compartment_unit = model.getVolumeUnits();
    } else if (compartment_list->get(species.getCompartment())->getSpatialDimensions() == 2) {
      compartment_unit = model.getAreaUnits(); 
    } else if (compartment_list->get(species.getCompartment())->getSpatialDimensions() == 1) {
      compartment_unit = model.getLengthUnits(); 
    }
    const LIBSBML_CPP_NAMESPACE::UnitDefinition* unit_definition = unit_definition_list->get(compartment_unit); 
    const LIBSBML_CPP_NAMESPACE::ListOfUnits* unit_list
    = unit_definition->getListOfUnits();
    unsigned int unitsSize = unit_list->size();
    double comp_unit = 1.0;
    for (unsigned int iu = 0u; iu < unitsSize; iu++) { 
      const LIBSBML_CPP_NAMESPACE::Unit* unit = unit_list->get(iu);
      comp_unit = comp_unit * pow(unit->getMultiplier()*pow(10,unit->getScale()),unit->getExponent());
    } 
    double avog_num = 6.02214e+23; 
    dynamic_cast<Species*>(m_p.get())->
      set_count(static_cast<species_cnt_t>(species.getInitialConcentration() * (6.02214E23 * 
      compartment_size) *comp_unit)); 
  }
}

template <typename G>
Vertex::Vertex(
  const LIBSBML_CPP_NAMESPACE::Model& model,
  const LIBSBML_CPP_NAMESPACE::Reaction& reaction,
  const G& g,
  const std::function<reaction_rate_t (const std::vector<reaction_rate_t>&)>&
  reaction_function_rate
  )
: m_type(_reaction_),
  m_typeid(static_cast<int>(_reaction_)),
  m_label(reaction.getIdAttribute()),
  m_pid(unassigned_partition),
  m_p(nullptr)
{
  using v_desc_t = typename boost::graph_traits<G>::vertex_descriptor;
  m_p = std::unique_ptr< Reaction<v_desc_t> >(new Reaction<v_desc_t>);
  dynamic_cast<Reaction<v_desc_t>*>(m_p.get())->set_rate_constant(1);

  using reaction_parameters_t = std::unordered_set<std::string>;
  typename reaction_parameters_t::const_iterator mpit, lpit;
  reaction_parameters_t mpset, lpset, pset;
  std::string formula;

  if (reaction.isSetKineticLaw()) {
    formula = SBML_formulaToString(reaction.getKineticLaw()->getMath());
  } else {
    WCS_THROW("The formula of the reaction " + reaction.getIdAttribute() + " should be set.");
  }
  //remove spaces from formula
  formula.erase(remove(formula.begin(), formula.end(), ' '), formula.end());

  std::string wholeformula("");
  //Add parameters
  const LIBSBML_CPP_NAMESPACE::ListOfParameters* parameter_list
    = model.getListOfParameters();
  unsigned int parametersSize = parameter_list->size(); 
  //Add local parameters
  if (model.getLevel() > 2) {
    const LIBSBML_CPP_NAMESPACE::ListOfLocalParameters* local_parameter_list
       = reaction.getKineticLaw()->getListOfLocalParameters();
    unsigned int num_local_parameters = reaction.getKineticLaw()->getNumLocalParameters();

    // Create an unordered_set for all local parameters of reaction
    for (unsigned int pi = 0u; pi < num_local_parameters; pi++) {
      const LIBSBML_CPP_NAMESPACE::LocalParameter* localparameter
      = local_parameter_list->get(pi);
      lpset.insert(localparameter->getIdAttribute());
    }
  } else {
    const LIBSBML_CPP_NAMESPACE::ListOfParameters* local_parameter_list
       = reaction.getKineticLaw()->getListOfParameters();
    unsigned int num_local_parameters = reaction.getKineticLaw()->getNumParameters();
    // Create an unordered_set for all local parameters of reaction
    for (unsigned int pi = 0u; pi < num_local_parameters; pi++) {
      const LIBSBML_CPP_NAMESPACE::Parameter* localparameter
      = local_parameter_list->get(pi);
      lpset.insert(localparameter->getIdAttribute());
    }
  }


  sbml_utils sbml_o;
  pset=sbml_o.get_reaction_parameters(model, reaction);

  // Create an unordered_set for all model parameters
  for (unsigned int pi = 0u; pi < parametersSize; pi++) {
    const LIBSBML_CPP_NAMESPACE::Parameter* parameter = parameter_list->get(pi);
    mpset.insert(parameter->getIdAttribute());
  }


  // Put initial assignments in map
  using  all_assignments_type
  = std::unordered_map<std::string, const LIBSBML_CPP_NAMESPACE::ASTNode *>;
  typename all_assignments_type::const_iterator allassigit;

  const LIBSBML_CPP_NAMESPACE::ListOfInitialAssignments* assignments_list
  = model.getListOfInitialAssignments();
  const unsigned int num_assignments = assignments_list->size();
  // A map for all assignments
  all_assignments_type all_assignments;
  //add initial assignments
  for (unsigned int ic = 0u; ic < num_assignments; ic++) {
    const LIBSBML_CPP_NAMESPACE::InitialAssignment& initial_assignment
    = *(assignments_list->get(ic));
    all_assignments.insert(std::make_pair(initial_assignment.getSymbol(),
                                          initial_assignment.getMath()));
  }
  //add rate rule assigments
  const LIBSBML_CPP_NAMESPACE::ListOfRules* rules_list = model.getListOfRules();
  const unsigned int rulesSize = rules_list->size();
  //add rate rules
  for (unsigned int ic = 0u; ic < rulesSize; ic++) {
    const LIBSBML_CPP_NAMESPACE::Rule& rule
    = *(rules_list->get(ic));
    all_assignments.insert(std::make_pair(rule.getVariable(),
                                          rule.getMath()));
  }

  //Add parameters
  for  (const std::string& x: pset) {
    std::string s_label = x;
    // model parameters
    mpit = mpset.find(s_label) ;
    if (mpit != mpset.end()) {
      std::stringstream ss;
      if (parameter_list->get(s_label)->isSetValue()){
        ss << parameter_list->get(s_label)->getValue();
        std::string parametervalue = ss.str();
        wholeformula = wholeformula + "var " + s_label
                     + " := " + parametervalue +  "; ";
      } else {
        allassigit = all_assignments.find(s_label);
        if (allassigit != all_assignments.end()) {
          ss << LIBSBML_CPP_NAMESPACE::SBML_formulaToString(allassigit->second);
          std::string parametervalue = ss.str();
          wholeformula = wholeformula + " " + s_label
                       + " := " + parametervalue +  "; ";
        }
      }
    }
    // reaction local parameters
    lpit = lpset.find(s_label) ;
    if (model.getLevel() > 2) {
      const LIBSBML_CPP_NAMESPACE::ListOfLocalParameters* local_parameter_list
      = reaction.getKineticLaw()->getListOfLocalParameters();
      if (lpit != lpset.end()) {
        std::stringstream ss;
        ss << local_parameter_list->get(s_label)->getValue();
        std::string parametervalue = ss.str();
        wholeformula = wholeformula + "var " + s_label
                    + " := " + parametervalue +  "; ";
      }
    } else {
      const LIBSBML_CPP_NAMESPACE::ListOfParameters* local_parameter_list
       = reaction.getKineticLaw()->getListOfParameters(); 
      if (lpit != lpset.end()) {
        std::stringstream ss;
        ss << local_parameter_list->get(s_label)->getValue();
        std::string parametervalue = ss.str();
        wholeformula = wholeformula + "var " + s_label
                    + " := " + parametervalue +  "; ";
      } 
    }  
  }

  //Add compartments
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
  if (reaction.isSetAnnotation()) {
    std::string annot= reaction.getAnnotationString() ;
    m_annotation = annot; 
  } else {
    m_annotation = ""; 
  }

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
    WCS_THROW("Attempted to dereference a wrong type of property pointer.");
  }
  return *ptr;
}

/**@}*/
} // end of namespace wcs

#endif // __WCS_REACTION_NETWORK_VERTEX_HPP__
