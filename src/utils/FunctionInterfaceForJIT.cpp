/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include "utils/FunctionInterfaceForJIT.hpp"

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

FunctionInterfaceForJIT::FunctionInterfaceForJIT()
{}

#if defined(WCS_HAS_SBML)
using ast_function_ptr = std::function<std::string (LIBSBML_CPP_NAMESPACE::ASTNode_t*)>;
using node_structure = std::tuple<ast_function_ptr, std::string, std::string, std::string>; 
using node_map = std::unordered_map<LIBSBML_CPP_NAMESPACE::ASTNodeType_t, node_structure>; 
node_map nodetypes = {
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_ABS,{nullptr,"std::abs(", "", ")"}},
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_ARCCOS,{nullptr,"std::acos(", "", ")"}},
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_ARCCOSH,{nullptr,"std::acosh(", "", ")"}},
        // {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_ARCCOT,{nullptr,"std::(", "", ")"}},  // arccotangent
        // {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_ARCCOTH,{nullptr,"std::(", "", ")"}}, // Hyperbolic arccotangent
        // {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_ARCCSC,{nullptr,"std::(", "", ")"}},  // arccosecant
        // {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_ARCCSCH,{nullptr,"std::(", "", ")"}}, // Hyperbolic arccosecant
        // {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_ARCSEC,{nullptr,"std::(", "", ")"}},  // arcsecant
        // {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_ARCSECH,{nullptr,"std::(", "", ")"}}, // Hyperbolic arcsecant
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_ARCSIN,{nullptr,"std::asin(", "", ")"}},
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_ARCSINH,{nullptr,"std::asinh(", "", ")"}},
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_ARCTAN,{nullptr,"std::atan(", "", ")"}},
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_ARCTANH,{nullptr,"std::atanh(", "", ")"}},
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_CEILING,{nullptr,"std::ceil(", "", ")"}},
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_COS,{nullptr,"std::cos(", "", ")"}},
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_COSH,{nullptr,"std::cosh(", "", ")"}},
        // {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_COT,{nullptr,"std::(", "", ")"}},  //cotangent
        // {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_COTH,{nullptr,"std::(", "", ")"}}, // Hyperbolic cotangent
        // {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_CSC,{nullptr,"std::(", "", ")"}},  // cosecant 
        // {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_CSCH,{nullptr,"std::(", "", ")"}}, // Hyperbolic cosecant
        // {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_DELAY,{nullptr,"std::(", "", ")"}},  //Delay
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_EXP,{nullptr,"std::exp(", "", ")"}},
        // {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_FACTORIAL,{nullptr,"std::(", "", ")"}}, // Factorial
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_FLOOR,{nullptr,"std::floor(", "", ")"}},
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_LN,{nullptr,"std::log(", "", ")"}},  // natural log
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_LOG,{nullptr,"", "", ""}}, // function pointer-------
        // {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_PIECEWISE,{nullptr,"std::(", "", ")"}}, // Piecewise
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_POWER,{nullptr,"std::pow(", ", ", ")"}}, // NOT SURE
        // {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_SEC,{nullptr,"std::(", "", ")"}}, // Secant
        // {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_SECH,{nullptr,"std::(", "", ")"}}, // Hyperbolic Secant
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_SIN,{nullptr,"std::sin(", "", ")"}},
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_SINH,{nullptr,"std::sinh(", "", ")"}},
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_TAN,{nullptr,"std::tan(", "", ")"}},
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_TANH,{nullptr,"std::tanh(", "", ")"}},

        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_MIN,{nullptr,"std::min({", ", ", "})"}},
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_MAX,{nullptr,"std::max({", ", ", "})"}},
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_QUOTIENT,{nullptr,"", "", ""}}, // function pointer div -------
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_RATE_OF,{nullptr,"", "", ""}}, // function pointer rate_of -------
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_REM,{nullptr,"std::remainder(", ", ", ")"}},
        
        {LIBSBML_CPP_NAMESPACE::AST_FUNCTION_ROOT,{nullptr,"", "", ""}}, //function pointer
    };
/**
 * Visits the given ASTNode as a function.  For this node only the
 * traversal is preorder.
 */
/* WCS addition: Add "std::" before min and "{...}" between parenthesis in order to support min */
void FunctionInterfaceForJIT::FormulaFormatter_visitFunction_wcs ( const ASTNode_t *parent,
                                 const ASTNode_t *node,
                                 strbuff_t  *sb )
{
  unsigned int numChildren = ASTNode_getNumChildren(node);
  unsigned int n;

  // int node_type = ASTNode_getType(node);
  LIBSBML_CPP_NAMESPACE::ASTNodeType_t node_type = ASTNode_getType(node);
  constexpr LIBSBML_CPP_NAMESPACE::ASTNodeType_t ast_func_min = LIBSBML_CPP_NAMESPACE::AST_FUNCTION_MIN;
  constexpr LIBSBML_CPP_NAMESPACE::ASTNodeType_t ast_func_max = LIBSBML_CPP_NAMESPACE::AST_FUNCTION_MAX;
  if (node_type == ast_func_max || node_type == ast_func_min) { 
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
    FunctionInterfaceForJIT::FormulaFormatter_visit_wcs( node, LIBSBML_CPP_NAMESPACE::ASTNode_getChild(node, 0), sb );
  }

  for (n = 1; n < numChildren; n++)
  {
    StringBuffer_appendChar(sb, ',');
    StringBuffer_appendChar(sb, ' ');
    FunctionInterfaceForJIT::FormulaFormatter_visit_wcs( node, LIBSBML_CPP_NAMESPACE::ASTNode_getChild(node, n), sb );
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
void FunctionInterfaceForJIT::FormulaFormatter_visitLog10_wcs ( const ASTNode_t *parent,
                              const ASTNode_t *node,
                              strbuff_t  *sb )
{
  StringBuffer_append(sb, "log10(");
  FunctionInterfaceForJIT::FormulaFormatter_visit_wcs(node, ASTNode_getChild(node, 1), sb);
  StringBuffer_appendChar(sb, ')');
}


/**
 * Visits the given ASTNode as the function "root(2, x)" and in doing so,
 * formats it as "sqrt(x)" (where x is any subexpression).
 */
/* WCS addition: Change function call from FormulaFormatter_visit to FormulaFormatter_visit_wcs  */
void FunctionInterfaceForJIT::FormulaFormatter_visitSqrt_wcs ( const ASTNode_t *parent,
                             const ASTNode_t *node,
                             strbuff_t  *sb )
{
  StringBuffer_append(sb, "sqrt(");
  FunctionInterfaceForJIT::FormulaFormatter_visit_wcs(node, ASTNode_getChild(node, 1), sb);
  StringBuffer_appendChar(sb, ')');
}


/**
 * Visits the given ASTNode as a unary minus.  For this node only the
 * traversal is preorder.
 */
/* WCS addition: Change function call from FormulaFormatter_visit to FormulaFormatter_visit_wcs  */
void FunctionInterfaceForJIT::FormulaFormatter_visitUMinus_wcs ( const ASTNode_t *parent,
                               const ASTNode_t *node,
                               strbuff_t  *sb )
{
  StringBuffer_appendChar(sb, '-');
  FunctionInterfaceForJIT::FormulaFormatter_visit_wcs ( node, ASTNode_getLeftChild(node), sb );
}


/**
 * Visits the given ASTNode and continues the inorder traversal.
 */
/* WCS addition: Change function call from FormulaFormatter_visit to FormulaFormatter_visit_wcs  */
void FunctionInterfaceForJIT::FormulaFormatter_visitOther_wcs ( const ASTNode_t *parent,
                              const ASTNode_t *node,
                              strbuff_t  *sb )
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
    FunctionInterfaceForJIT::FormulaFormatter_visit_wcs( node, ASTNode_getChild(node, 0), sb );
    StringBuffer_appendChar(sb, ')');
  }

  else {
    FunctionInterfaceForJIT::FormulaFormatter_visit_wcs( node, ASTNode_getChild(node, 0), sb );

    for (n = 1; n < numChildren; n++)
    {
      FormulaFormatter_format(sb, node);
      FunctionInterfaceForJIT::FormulaFormatter_visit_wcs( node, ASTNode_getChild(node, n), sb );
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
void FunctionInterfaceForJIT::FormulaFormatter_visit_wcs ( const ASTNode_t *parent,
                         const ASTNode_t *node,
                         strbuff_t  *sb )
{
  if (ASTNode_isLog10(node))
  {
    StringBuffer_append(sb, "log10(");
    FunctionInterfaceForJIT::FormulaFormatter_visit_wcs(node, ASTNode_getChild(node, 1), sb);
    StringBuffer_appendChar(sb, ')');
    // FormulaFormatter_visitLog10_wcs(parent, node, sb);
  }
  else if (ASTNode_isSqrt(node))
  {
    StringBuffer_append(sb, "sqrt(");
    FunctionInterfaceForJIT::FormulaFormatter_visit_wcs(node, ASTNode_getChild(node, 1), sb);
    StringBuffer_appendChar(sb, ')');
    // FormulaFormatter_visitSqrt_wcs(parent, node, sb);
  }
  else if (FormulaFormatter_isFunction(node))
  {
    node_map::iterator nit;
    LIBSBML_CPP_NAMESPACE::ASTNodeType_t node_type = ASTNode_getType(node);
    nit = nodetypes.find(node_type);
    if (nit != nodetypes.cend()) {
      node_structure node_st = nit -> second;
      ast_function_ptr node_f_ptr = std::get<0>(node_st); 
      if (node_f_ptr == nullptr) {
        unsigned int numChildren = ASTNode_getNumChildren(node);
        unsigned int n;
        std::string fist_symbol = std::get<1>(node_st);
        std::string middle_symbol = std::get<2>(node_st);
        std::string end_symbol = std::get<3>(node_st);
        StringBuffer_append(sb, fist_symbol.c_str());
        // FormulaFormatter_format(sb, node);
        if (numChildren > 0) {
          FunctionInterfaceForJIT::FormulaFormatter_visit_wcs( node, LIBSBML_CPP_NAMESPACE::ASTNode_getChild(node, 0), sb );
        }
        for (n = 1; n < numChildren; n++) {
          StringBuffer_append(sb, middle_symbol.c_str()); 
          FunctionInterfaceForJIT::FormulaFormatter_visit_wcs( node, LIBSBML_CPP_NAMESPACE::ASTNode_getChild(node, n), sb );
        }
        StringBuffer_append(sb, end_symbol.c_str()); 
      } else { // if there is a function pointer

      }   
     } else { // user functions
        unsigned int numChildren = ASTNode_getNumChildren(node);
        unsigned int n;
        StringBuffer_appendChar(sb, '(');
        StringBuffer_append(sb, ASTNode_getName(node));

        if (numChildren > 0)
        {
          FunctionInterfaceForJIT::FormulaFormatter_visit_wcs( node, LIBSBML_CPP_NAMESPACE::ASTNode_getChild(node, 0), sb );
        }

        for (n = 1; n < numChildren; n++)
        {
          StringBuffer_appendChar(sb, ',');
          StringBuffer_appendChar(sb, ' ');
          FunctionInterfaceForJIT::FormulaFormatter_visit_wcs( node, LIBSBML_CPP_NAMESPACE::ASTNode_getChild(node, n), sb );
        }
        StringBuffer_appendChar(sb, ')');
    }
    // FormulaFormatter_visitFunction_wcs(parent, node, sb);
  }
  else if (ASTNode_hasTypeAndNumChildren(node, LIBSBML_CPP_NAMESPACE::AST_MINUS, 1))
  {
    StringBuffer_appendChar(sb, '-');
    FunctionInterfaceForJIT::FormulaFormatter_visit_wcs ( node, ASTNode_getLeftChild(node), sb );
    // FormulaFormatter_visitUMinus_wcs(parent, node, sb);
  }
  else if (ASTNode_hasTypeAndNumChildren(node, LIBSBML_CPP_NAMESPACE::AST_PLUS, 1) || ASTNode_hasTypeAndNumChildren(node, LIBSBML_CPP_NAMESPACE::AST_TIMES, 1))
  {
    FunctionInterfaceForJIT::FormulaFormatter_visit_wcs (node, ASTNode_getChild(node, 0), sb);
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
std::string FunctionInterfaceForJIT::SBML_formulaToString_wcs (const ASTNode_t *tree)
{
  if (tree == nullptr) {
    return nullptr;
  }
  
  LIBSBML_CPP_NAMESPACE::StringBuffer_t *buf = LIBSBML_CPP_NAMESPACE::StringBuffer_create(1024);

  FormulaFormatter_visit_wcs(NULL, tree, buf);
  std::string str(LIBSBML_CPP_NAMESPACE::StringBuffer_getBuffer(buf));
  safe_free(buf);
  
  return str;
}
#endif // defined(WCS_HAS_SBML)

FunctionInterfaceForJIT::~FunctionInterfaceForJIT()
{}

/**@}*/
} // end of namespace wcs
