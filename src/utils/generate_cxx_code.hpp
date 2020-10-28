/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef  __WCS_UTILS_GENERATE_CXX_CODE__
#define  __WCS_UTILS_GENERATE_CXX_CODE__

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#include <unordered_map>
#include "sbml_utils.hpp"
#include <cstdio>

#if defined(WCS_HAS_SBML)
#include <sbml/SBMLTypes.h>
#include <sbml/common/extern.h>

LIBSBML_CPP_NAMESPACE_USE

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

class generate_cxx_code {
 public:
  using map_symbol_to_ast_node_t
    = std::unordered_map <std::string, const LIBSBML_CPP_NAMESPACE::ASTNode *>;

  static const std::string generate_code(const LIBSBML_CPP_NAMESPACE::Model& model);
  static const std::string compile_code(const std::string generated_filename);

  template <typename TTT>
  class basetype_to_string {
   public:
    static const char* value;
  };

 private:
  static void get_dependencies(
    const LIBSBML_CPP_NAMESPACE::ASTNode& math,
    std::vector<std::string> & math_elements,
    std::unordered_set<std::string> & good,
    map_symbol_to_ast_node_t & constant_init_assig,
    map_symbol_to_ast_node_t & variables_assig_rul,
    map_symbol_to_ast_node_t & model_reactions_s);

  static std::vector<std::string> get_all_dependencies(
    const LIBSBML_CPP_NAMESPACE::ASTNode& formula,
    std::unordered_set<std::string> & good,
    map_symbol_to_ast_node_t & constant_init_assig,
    map_symbol_to_ast_node_t & variables_assig_rul,
    map_symbol_to_ast_node_t & model_reactions_s);

};

/**@}*/
} // end of namespace wcs
#endif // defined(WCS_HAS_SBML)
#endif //  __WCS_UTILS_GENERATE_CXX_CODE__