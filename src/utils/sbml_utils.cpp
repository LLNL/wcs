/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include "utils/sbml_utils.hpp"

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif


#if defined(WCS_HAS_SBML)
#include <sbml/SBMLTypes.h>
#include <sbml/common/extern.h>
#endif // defined(WCS_HAS_SBML)

#include <iostream>
#include <string>
#include <unordered_set>
#include <utils/exception.hpp>

namespace wcs {
/** \addtogroup wcs_utils
 *  *  @{ */

sbml_utils::sbml_utils()
{}

#if defined(WCS_HAS_SBML)

std::unordered_set<std::string> sbml_utils::find_undeclared_species_in_reaction_formula(const LIBSBML_CPP_NAMESPACE
::Model& model, const LIBSBML_CPP_NAMESPACE::Reaction& reaction)
{
  using model_parameters = std::unordered_set<std::string>;
  model_parameters pset;

  using reaction_reactants = std::unordered_set<std::string>;
  reaction_reactants rset;

  reaction_reactants undeclared_reactants;

  using reaction_modifiers = std::unordered_set<std::string>;
  reaction_modifiers mset;


  using model_compartments = std::unordered_set<std::string>;
  model_compartments cset;

  using reaction_symbols = std::unordered_set<std::string>;
  reaction_symbols symbol_set;


  const LIBSBML_CPP_NAMESPACE::ListOfSpecies* specieslist = model.getListOfSpecies();

  if (specieslist == nullptr) {
    WCS_THROW("The species list of the SBML file is empty.");
  }

  // create an unordered_set for model parameters
  const LIBSBML_CPP_NAMESPACE::ListOfParameters* parameterslist = model.getListOfParameters();
  unsigned int parametersSize = parameterslist->size();

  if (parameterslist == nullptr) {
    WCS_THROW("The parameter list of the SBML file is empty.");
  }
  for (unsigned int pi = 0u; pi < parametersSize; pi++) {
    pset.insert(parameterslist->get(pi)->getIdAttribute());
  }

  // create an unordered set for reaction_reactants
  unsigned int reactantsSize = reaction.getNumReactants();
  for (unsigned int ri = 0u; ri < reactantsSize; ri++) {
      const auto &reactant = *(reaction.getReactant(ri));
      rset.insert(specieslist->get(reactant.getSpecies())->getIdAttribute());
  }

  // create an unordered set for reaction_modifiers
  unsigned int modifiersSize = reaction.getNumModifiers();
  for (unsigned int mi = 0u; mi < modifiersSize; mi++) {
      const auto &modifier = *(reaction.getModifier(mi));
      mset.insert(specieslist->get(modifier.getSpecies())->getIdAttribute());
  }

  // create an unordered_set for model compartments
  const LIBSBML_CPP_NAMESPACE::ListOfCompartments* compartmentslist = model.getListOfCompartments();
  unsigned int compartmentsSize = compartmentslist->size();
  if (compartmentslist == nullptr) {
    WCS_THROW("The compartment list of the SBML file is empty.");
  }
  for (unsigned int ci = 0u; ci < compartmentsSize; ci++) {
    cset.insert(compartmentslist->get(ci)->getIdAttribute());
  }

  const LIBSBML_CPP_NAMESPACE::ASTNode* formula = reaction.getKineticLaw()->getMath();
  if (formula == nullptr) {
    WCS_THROW("The " + reaction.getIdAttribute() + " reaction of the SBML file does not have a rate formula.");
  }
  symbol_set = get_symbol_table_of_formula(*formula);
  for  (const std::string& x: symbol_set) {
    if (pset.find(x) == pset.cend() && rset.find(x) == rset.cend() && cset.find(x) == cset.cend() && mset.find(x) == mset.cend()) {
      undeclared_reactants.insert(x);
    }
  }
  /**for  (const std::string& x: undeclared_reactants) {
   std::cout << " " << x;
  }
  std::cout <<  std::endl;*/
  return undeclared_reactants;
}


std::unordered_set<std::string> sbml_utils::get_reaction_parameters(const LIBSBML_CPP_NAMESPACE::Model& model,
 const LIBSBML_CPP_NAMESPACE::Reaction& reaction)
{
  using reaction_parameters = std::unordered_set<std::string>;
  reaction_parameters pset;

  using reaction_reactants= std::unordered_set<std::string>;
  typename reaction_reactants::const_iterator rit;
  reaction_reactants rset;


  using model_compartments= std::unordered_set<std::string>;
  typename model_compartments::const_iterator cit;
  model_compartments cset;

  using reaction_symbols= std::unordered_set<std::string>;
  reaction_symbols symbol_set;


  const LIBSBML_CPP_NAMESPACE::ListOfSpecies* specieslist = model.getListOfSpecies();


  // create an unordered set for reaction_reactants
  unsigned int reactantsSize = reaction.getNumReactants();
  for (unsigned int ri = 0u; ri < reactantsSize; ri++) {
      const auto &reactant = *(reaction.getReactant(ri));
      rset.insert(specieslist->get(reactant.getSpecies())->getIdAttribute());
  }

  // create an unordered_set for model compartments
  const LIBSBML_CPP_NAMESPACE::ListOfCompartments* compartmentslist = model.getListOfCompartments();
  unsigned int compartmentsSize = compartmentslist->size();
  for (unsigned int ci = 0u; ci < compartmentsSize; ci++) {
    cset.insert(compartmentslist->get(ci)->getIdAttribute());
  }

  const LIBSBML_CPP_NAMESPACE::ASTNode* formula = reaction.getKineticLaw()->getMath();
  symbol_set = get_symbol_table_of_formula(*formula);
  for  (const std::string& x: symbol_set) {
    rit = rset.find(x);
    cit = cset.find(x);
    if (rit == rset.end() && cit == cset.end()) {
      pset.insert(x);
    }
  }
  return pset;
}



std::unordered_set<std::string>  sbml_utils::get_symbol_table_of_formula(const LIBSBML_CPP_NAMESPACE::ASTNode& formula)
{
  std::unordered_set<std::string> symbol_table;
  read_ast_node(formula, symbol_table);
  /**for  (const std::string& x: symbol_table) {
   std::cout << " " << x;
  }
  std::cout <<  std::endl;*/
  return symbol_table;
}

void sbml_utils::read_ast_node(const LIBSBML_CPP_NAMESPACE::ASTNode&  math, std::unordered_set<std::string> & math_elements )
{
  if (math.getNumChildren() >= 1u) {
    for (unsigned int n = 0u; n < math.getNumChildren(); ++n) {
        const LIBSBML_CPP_NAMESPACE::ASTNode*  math_c = math.getChild(n);
        read_ast_node(*math_c, math_elements);
    }
  } else {
    if (!math.isNumber()) {
      char* formula = SBML_formulaToString(&math);
      std::unordered_set<std::string>::const_iterator mit = math_elements.find (formula);
      if (mit == math_elements.cend()) {
        math_elements.insert(formula);
      }
    }
  }
}
#endif // defined(WCS_HAS_SBML)

sbml_utils::~sbml_utils()
{}

/**@}*/
} // end of namespace wcs