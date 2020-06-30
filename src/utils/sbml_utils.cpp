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

namespace wcs {
/** \addtogroup wcs_utils
 *  *  @{ */

sbml_utils::sbml_utils()
: m_type(_undefined_)
{
}

#if defined(WCS_HAS_SBML)

std::unordered_set<std::string> sbml_utils::find_undeclared_reactant_species(const libsbml::Model& model, const libsbml::Reaction& reaction)
{
  using model_parameters = std::unordered_set<std::string>;
  typename model_parameters::const_iterator pit;
  model_parameters pset;

  using reaction_reactants= std::unordered_set<std::string>;
  typename reaction_reactants::const_iterator rit;
  reaction_reactants rset;

  reaction_reactants undeclared_reactants;
  
  using model_compartments= std::unordered_set<std::string>;
  typename model_compartments::const_iterator cit;
  model_compartments cset;

  using reaction_symbols= std::unordered_set<std::string>;
  ///typename reaction_symbols::const_iterator sit;
  reaction_symbols symbol_set;
  

  const libsbml::ListOfSpecies* specieslist = model.getListOfSpecies();
  
  // create an unordered_set for model parameters
  const libsbml::ListOfParameters* parameterslist = model.getListOfParameters();
  unsigned int parametersSize = parameterslist->size();
  for (unsigned int pi = 0u; pi < parametersSize; pi++) {
    pset.insert(parameterslist->get(pi)->getIdAttribute());
  }

  // create an unordered set for reaction_reactants
  unsigned int reactantsSize = reaction.getNumReactants();
  for (unsigned int ri = 0u; ri < reactantsSize; ri++) {
      const auto &reactant = *(reaction.getReactant(ri));
      rset.insert(specieslist->get(reactant.getSpecies())->getIdAttribute());
  }

  // create an unordered_set for model compartments
  const libsbml::ListOfCompartments* compartmentslist = model.getListOfCompartments();
  unsigned int compartmentsSize = compartmentslist->size();
  for (unsigned int ci = 0u; ci < compartmentsSize; ci++) {
    cset.insert(compartmentslist->get(ci)->getIdAttribute());
  }

  const libsbml::ASTNode* formula = reaction.getKineticLaw()->getMath();
  symbol_set = get_symbol_table_of_formula(*formula);
  for  (const std::string& x: symbol_set) {
    pit = pset.find(x);
    rit = rset.find(x);
    cit = cset.find(x);
    if (pit == pset.end() && rit == rset.end() && cit == cset.end()) {
      undeclared_reactants.insert(x);
    } 
  }   
  /**for  (const std::string& x: undeclared_reactants) {
   std::cout << " " << x;
  }
  std::cout <<  std::endl;*/
  return undeclared_reactants;
}


std::unordered_set<std::string>  sbml_utils::get_symbol_table_of_formula(const libsbml::ASTNode& formula) 
{
  std::unordered_set<std::string> symbol_table;
  read_ast_node(formula, symbol_table);
  return symbol_table;
}

void sbml_utils::read_ast_node(const libsbml::ASTNode&  math, std::unordered_set<std::string> & math_elements )
{
  if (math.getNumChildren() > 1) {
    for (unsigned int n = 0; n < math.getNumChildren(); ++n) {
        const libsbml::ASTNode*  math_c = math.getChild(n);
        read_ast_node(*math_c, math_elements);
    }
  } else {
    if (!math.isNumber()) {
      char* formula = SBML_formulaToString(&math);
      std::unordered_set<std::string>::const_iterator mit = math_elements.find (formula);
      if (mit == math_elements.end()) {      
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