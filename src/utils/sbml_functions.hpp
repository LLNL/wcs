/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef  __WCS_UTILS_SBML_FUNCTIONS__
#define  __WCS_UTILS_SBML_FUNCTIONS__

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
/** \addtogroup wcs_functions
 *  @{ */

class sbml_functions {
 public:
  sbml_functions();
  ~sbml_functions();

  #if defined(WCS_HAS_SBML)
  static void
  FormulaFormatter_visit_wcs ( const LIBSBML_CPP_NAMESPACE::ASTNode_t *parent,
                         const LIBSBML_CPP_NAMESPACE::ASTNode_t *node,
                         LIBSBML_CPP_NAMESPACE::StringBuffer_t  *sb );

  
  static char *
  SBML_formulaToString_wcs (const LIBSBML_CPP_NAMESPACE::ASTNode_t *tree);

  static void
  FormulaFormatter_visitFunction_wcs ( const LIBSBML_CPP_NAMESPACE::ASTNode_t *parent,
                                  const LIBSBML_CPP_NAMESPACE::ASTNode_t *node,
                                  LIBSBML_CPP_NAMESPACE::StringBuffer_t  *sb );

  static void
  FormulaFormatter_visitLog10_wcs ( const LIBSBML_CPP_NAMESPACE::ASTNode_t *parent,
                              const LIBSBML_CPP_NAMESPACE::ASTNode_t *node,
                              LIBSBML_CPP_NAMESPACE::StringBuffer_t  *sb );

  static void
  FormulaFormatter_visitSqrt_wcs ( const LIBSBML_CPP_NAMESPACE::ASTNode_t *parent,
                             const LIBSBML_CPP_NAMESPACE::ASTNode_t *node,
                             LIBSBML_CPP_NAMESPACE::StringBuffer_t  *sb );

  static void
  FormulaFormatter_visitUMinus_wcs ( const LIBSBML_CPP_NAMESPACE::ASTNode_t *parent,
                               const LIBSBML_CPP_NAMESPACE::ASTNode_t *node,
                               LIBSBML_CPP_NAMESPACE::StringBuffer_t  *sb );

  static void
  FormulaFormatter_visitOther_wcs ( const LIBSBML_CPP_NAMESPACE::ASTNode_t *parent,
                              const LIBSBML_CPP_NAMESPACE::ASTNode_t *node,
                              LIBSBML_CPP_NAMESPACE::StringBuffer_t  *sb );

  #endif // defined(WCS_HAS_SBML)

};

/**@}*/
} // end of namespace wcs
#endif //  __WCS_UTILS_SBML_FUNCTIONS__
