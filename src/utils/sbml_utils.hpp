/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef  __WCS_UTILS_SBML_UTILS__
#define  __WCS_UTILS_SBML_UTILS__

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#include <string>
#include <unordered_set>

#if defined(WCS_HAS_SBML)
#include <sbml/SBMLTypes.h>
#include <sbml/common/extern.h>
#endif // defined(WCS_HAS_SBML)


namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

class sbml_utils {
 public:
  sbml_utils();
  ~sbml_utils();

  #if defined(WCS_HAS_SBML)
  static std::unordered_set<std::string>
  find_undeclared_species_in_reaction_formula(
    const LIBSBML_CPP_NAMESPACE::Model& model,
    const LIBSBML_CPP_NAMESPACE::Reaction& reaction);

  static std::unordered_set<std::string>
  get_reaction_parameters(
    const LIBSBML_CPP_NAMESPACE::Model& model,
    const LIBSBML_CPP_NAMESPACE::Reaction& reaction);

  static std::unordered_set<std::string>
  get_symbol_table_of_formula(
    const LIBSBML_CPP_NAMESPACE::ASTNode& formula);

  static void
  read_ast_node(
    const LIBSBML_CPP_NAMESPACE::ASTNode& math,
    std::unordered_set<std::string>& math_elements);

  #endif // defined(WCS_HAS_SBML)

};

/**@}*/
} // end of namespace wcs
#endif //  __WCS_UTILS_SBML_UTILS__
