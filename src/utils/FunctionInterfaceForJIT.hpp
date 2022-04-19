/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef  __WCS_UTILS_FUNCTIONINTERFACEFORJIT__
#define  __WCS_UTILS_FUNCTIONINTERFACEFORJIT__

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

#if defined(WCS_HAS_SBML)
class FunctionInterfaceForJIT {
  public:
    using ASTNode_t = LIBSBML_CPP_NAMESPACE::ASTNode_t;
    using strbuff_t = LIBSBML_CPP_NAMESPACE::StringBuffer_t;
  public:
    FunctionInterfaceForJIT();
    ~FunctionInterfaceForJIT();

    static std::string SBML_formulaToString_wcs (const ASTNode_t *tree);
   
  private:
    static void FormulaFormatter_visit_wcs ( const ASTNode_t*parent,
                          const ASTNode_t *node,
                          strbuff_t  *sb );

    static void FormulaFormatter_visitFunction_wcs ( const ASTNode_t *parent,
                                    const ASTNode_t *node,
                                    strbuff_t  *sb );

    static void FormulaFormatter_visitLog10_wcs ( const ASTNode_t*parent,
                                const ASTNode_t *node,
                                strbuff_t *sb );

    static void FormulaFormatter_visitSqrt_wcs ( const ASTNode_t *parent,
                              const ASTNode_t *node,
                              strbuff_t  *sb );

    static void FormulaFormatter_visitUMinus_wcs ( const ASTNode_t *parent,
                                const ASTNode_t *node,
                                strbuff_t  *sb );

    static void FormulaFormatter_visitOther_wcs ( const ASTNode_t *parent,
                                const ASTNode_t *node,
                                strbuff_t *sb );

};
#endif // defined(WCS_HAS_SBML)

/**@}*/
} // end of namespace wcs
#endif //  __WCS_UTILS_FUNCTIONINTERFACEFORJIT__
