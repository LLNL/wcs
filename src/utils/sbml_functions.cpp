/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include "utils/sbml_functions.hpp"

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
// #include <unordered_set>
// #include <utils/exception.hpp>

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

sbml_functions::sbml_functions()
{}

#if defined(WCS_HAS_SBML)

/**
 * Visits the given ASTNode as a function.  For this node only the
 * traversal is preorder.
 */
/* WCS addition: Add "std::" before min and "{...}" between parenthesis in order to support min */
void
sbml_functions::FormulaFormatter_visitFunction_wcs ( const LIBSBML_CPP_NAMESPACE::ASTNode_t *parent,
                                 const LIBSBML_CPP_NAMESPACE::ASTNode_t *node,
                                 LIBSBML_CPP_NAMESPACE::StringBuffer_t  *sb )
{
  unsigned int numChildren = ASTNode_getNumChildren(node);
  unsigned int n;

  int node_type = ASTNode_getType(node);
  //321 for min, 320 fo rmax 
  // if (node_name == "min" || node_name == "max") {
  if (node_type == 321 || node_type == 320) { 
    StringBuffer_append(sb, "std::");
  }
  FormulaFormatter_format(sb, node);

  // support min in sbml 
  StringBuffer_appendChar(sb, '(');
  if (node_type == 321 || node_type == 320) { 
    StringBuffer_appendChar(sb, '{');
  }
  
  if (numChildren > 0)
  {
    sbml_functions::FormulaFormatter_visit_wcs( node, LIBSBML_CPP_NAMESPACE::ASTNode_getChild(node, 0), sb );
  }

  for (n = 1; n < numChildren; n++)
  {
    StringBuffer_appendChar(sb, ',');
    StringBuffer_appendChar(sb, ' ');
    sbml_functions::FormulaFormatter_visit_wcs( node, LIBSBML_CPP_NAMESPACE::ASTNode_getChild(node, n), sb );
  }

  if (node_type == 321 || node_type == 320) { 
    StringBuffer_appendChar(sb, '}');
  }
  StringBuffer_appendChar(sb, ')');
}


/**
 * Visits the given ASTNode as the function "log(10, x)" and in doing so,
 * formats it as "log10(x)" (where x is any subexpression).
 */
/* WCS addition: Change function call from FormulaFormatter_visit to FormulaFormatter_visit_wcs  */
void
sbml_functions::FormulaFormatter_visitLog10_wcs ( const LIBSBML_CPP_NAMESPACE::ASTNode_t *parent,
                              const LIBSBML_CPP_NAMESPACE::ASTNode_t *node,
                              LIBSBML_CPP_NAMESPACE::StringBuffer_t  *sb )
{
  StringBuffer_append(sb, "log10(");
  sbml_functions::FormulaFormatter_visit_wcs(node, ASTNode_getChild(node, 1), sb);
  StringBuffer_appendChar(sb, ')');
}


/**
 * Visits the given ASTNode as the function "root(2, x)" and in doing so,
 * formats it as "sqrt(x)" (where x is any subexpression).
 */
/* WCS addition: Change function call from FormulaFormatter_visit to FormulaFormatter_visit_wcs  */
void
sbml_functions::FormulaFormatter_visitSqrt_wcs ( const LIBSBML_CPP_NAMESPACE::ASTNode_t *parent,
                             const LIBSBML_CPP_NAMESPACE::ASTNode_t *node,
                             LIBSBML_CPP_NAMESPACE::StringBuffer_t  *sb )
{
  StringBuffer_append(sb, "sqrt(");
  sbml_functions::FormulaFormatter_visit_wcs(node, ASTNode_getChild(node, 1), sb);
  StringBuffer_appendChar(sb, ')');
}


/**
 * Visits the given ASTNode as a unary minus.  For this node only the
 * traversal is preorder.
 */
/* WCS addition: Change function call from FormulaFormatter_visit to FormulaFormatter_visit_wcs  */
void
sbml_functions::FormulaFormatter_visitUMinus_wcs ( const LIBSBML_CPP_NAMESPACE::ASTNode_t *parent,
                               const LIBSBML_CPP_NAMESPACE::ASTNode_t *node,
                               LIBSBML_CPP_NAMESPACE::StringBuffer_t  *sb )
{
  StringBuffer_appendChar(sb, '-');
  sbml_functions::FormulaFormatter_visit_wcs ( node, ASTNode_getLeftChild(node), sb );
}


/**
 * Visits the given ASTNode and continues the inorder traversal.
 */
/* WCS addition: Change function call from FormulaFormatter_visit to FormulaFormatter_visit_wcs  */
void
sbml_functions::FormulaFormatter_visitOther_wcs ( const LIBSBML_CPP_NAMESPACE::ASTNode_t *parent,
                              const LIBSBML_CPP_NAMESPACE::ASTNode_t *node,
                              LIBSBML_CPP_NAMESPACE::StringBuffer_t  *sb )
{
  unsigned int numChildren = ASTNode_getNumChildren(node);
  int group       = FormulaFormatter_isGrouped(parent, node);
  unsigned int n;


  if (group)
  {
    StringBuffer_appendChar(sb, '(');
  }

  if (numChildren == 0) {
    FormulaFormatter_format(sb, node);
  }

  else if (numChildren == 1)
  {
    //I believe this would only be called for invalid ASTNode setups,
    // but this could in theory occur.  This is the safest 
    // behavior I can think of.
    FormulaFormatter_format(sb, node);
    StringBuffer_appendChar(sb, '(');
    sbml_functions::FormulaFormatter_visit_wcs( node, ASTNode_getChild(node, 0), sb );
    StringBuffer_appendChar(sb, ')');
  }

  else {
    sbml_functions::FormulaFormatter_visit_wcs( node, ASTNode_getChild(node, 0), sb );

    for (n = 1; n < numChildren; n++)
    {
      FormulaFormatter_format(sb, node);
      sbml_functions::FormulaFormatter_visit_wcs( node, ASTNode_getChild(node, n), sb );
    }
  }

  if (group)
  {
    StringBuffer_appendChar(sb, ')');
  }
}



/**
 * Visits the given ASTNode node.  This function is really just a
 * dispatcher to either SBML_formulaToString_visitFunction() or
 * SBML_formulaToString_visitOther().
 */
/* WCS addition: Change function calls to wcs modified fuctions */
void
sbml_functions::FormulaFormatter_visit_wcs ( const LIBSBML_CPP_NAMESPACE::ASTNode_t *parent,
                         const LIBSBML_CPP_NAMESPACE::ASTNode_t *node,
                         LIBSBML_CPP_NAMESPACE::StringBuffer_t  *sb )
{
  if (ASTNode_isLog10(node))
  {
    FormulaFormatter_visitLog10_wcs(parent, node, sb);
  }
  else if (ASTNode_isSqrt(node))
  {
    FormulaFormatter_visitSqrt_wcs(parent, node, sb);
  }
  else if (FormulaFormatter_isFunction(node))
  {
    FormulaFormatter_visitFunction_wcs(parent, node, sb);
  }
  else if (ASTNode_hasTypeAndNumChildren(node, LIBSBML_CPP_NAMESPACE::AST_MINUS, 1))
  {
    FormulaFormatter_visitUMinus_wcs(parent, node, sb);
  }
  else if (ASTNode_hasTypeAndNumChildren(node, LIBSBML_CPP_NAMESPACE::AST_PLUS, 1) || ASTNode_hasTypeAndNumChildren(node, LIBSBML_CPP_NAMESPACE::AST_TIMES, 1))
  {
    libsbml::FormulaFormatter_visit(node, ASTNode_getChild(node, 0), sb);
  }
  else if (ASTNode_hasTypeAndNumChildren(node, LIBSBML_CPP_NAMESPACE::AST_PLUS, 0))
  {
    libsbml::StringBuffer_appendInt(sb, 0);
  }
  else if (ASTNode_hasTypeAndNumChildren(node, LIBSBML_CPP_NAMESPACE::AST_TIMES, 0))
  {
    libsbml::StringBuffer_appendInt(sb, 1);
  }
  else
  {
    FormulaFormatter_visitOther_wcs(parent, node, sb);
  }
}

/* WCS addition: Change function call from FormulaFormatter_visit to FormulaFormatter_visit_wcs  */
char *
sbml_functions::SBML_formulaToString_wcs (const LIBSBML_CPP_NAMESPACE::ASTNode_t *tree)
{
  char *s;

  if (tree == NULL)
  {
    s = NULL;
  }
  else
  {
    LIBSBML_CPP_NAMESPACE::StringBuffer_t *sb = LIBSBML_CPP_NAMESPACE::StringBuffer_create(128);

    FormulaFormatter_visit_wcs(NULL, tree, sb);
    s = LIBSBML_CPP_NAMESPACE::StringBuffer_getBuffer(sb);
    safe_free(sb);
  }
  return s;
}
#endif // defined(WCS_HAS_SBML)

sbml_functions::~sbml_functions()
{}

/**@}*/
} // end of namespace wcs
