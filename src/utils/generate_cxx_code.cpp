/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#include "utils/generate_cxx_code.hpp"
#include "utils/exception.hpp"
#include "utils/file.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h> //system
#include "wcs_types.hpp"
#include <set>
#include <regex>


#include <string>
#include <filesystem>

//#include <stdio.h>

#if defined(WCS_HAS_SBML)

LIBSBML_CPP_NAMESPACE_USE

namespace wcs {

using map_symbol_to_ast_node_t
  = std::unordered_map <std::string, const LIBSBML_CPP_NAMESPACE::ASTNode *>;

template<>
const char* generate_cxx_code::basetype_to_string<double>::value = "double";
template<>
const char* generate_cxx_code::basetype_to_string<float>::value = "float";
typedef reaction_rate_t ( * rate_function_pointer)(const std::vector<reaction_rate_t>&);

void
generate_cxx_code::get_dependencies(
  const LIBSBML_CPP_NAMESPACE::ASTNode& math,
  std::vector<std::string> & math_elements,
  std::unordered_set<std::string> & good,
  map_symbol_to_ast_node_t & constant_init_assig,
  map_symbol_to_ast_node_t & variables_assig_rul,
  map_symbol_to_ast_node_t & model_reactions_s)
{
  if (math.getNumChildren() >= 1u) {
    for (unsigned int n = 0u; n < math.getNumChildren(); ++n) {
      const LIBSBML_CPP_NAMESPACE::ASTNode*  math_c = math.getChild(n);
      get_dependencies(*math_c,
                       math_elements,
                       good,
                       constant_init_assig,
                       variables_assig_rul,
                       model_reactions_s);
    }
  } else {
    if (!math.isNumber()) {
      const char* formula = SBML_formulaToString(&math);
      map_symbol_to_ast_node_t::const_iterator inasit, asruit, mreactit;
      inasit = constant_init_assig.find(formula);
      asruit = variables_assig_rul.find(formula);
      mreactit = model_reactions_s.find(formula);
      std::unordered_set<std::string>::const_iterator git = good.find(formula);
      std::unordered_set<std::string> no_duplicates_vector;
      //if (constant_init_assig.empty() && variables_assig_rul.empty()){

      //} else {
      if (git == good.cend()) {
        if (inasit != constant_init_assig.cend()) {
          math_elements.insert(math_elements.begin(), formula);
          get_dependencies(*inasit->second,
                           math_elements,
                           good,
                           constant_init_assig,
                           variables_assig_rul,
                           model_reactions_s);

        } else if (asruit != variables_assig_rul.cend()) {
          math_elements.insert(math_elements.begin(), formula);
          get_dependencies(*asruit->second,
                           math_elements,
                           good,
                           constant_init_assig,
                           variables_assig_rul,
                           model_reactions_s);

        } else if (mreactit != model_reactions_s.cend()){
          math_elements.insert(math_elements.begin(), formula);
          get_dependencies(*mreactit->second,
                           math_elements,
                           good,
                           constant_init_assig,
                           variables_assig_rul,
                           model_reactions_s);

        } else {
          math_elements.insert(math_elements.begin(), formula);
        }
      //}
      }
    }
  }
}


std::vector<std::string>
generate_cxx_code::get_all_dependencies(
  const LIBSBML_CPP_NAMESPACE::ASTNode& formula,
  std::unordered_set<std::string> & good,
  map_symbol_to_ast_node_t & constant_init_assig,
  map_symbol_to_ast_node_t & variables_assig_rul,
  map_symbol_to_ast_node_t & model_reactions_s)
{

  std::vector<std::string> dependencies, dependencies_no_dupl;
  std::unordered_set<std::string> dependencies_set;
  std::unordered_set<std::string>::const_iterator ndit;

  get_dependencies(formula,
                   dependencies,
                   good,
                   constant_init_assig,
                   variables_assig_rul,
                   model_reactions_s);

  std::vector<std::string>::const_iterator it;
  //remove duplicates
  for  (it = dependencies.cbegin(); it < dependencies.cend(); it++) {
    ndit = dependencies_set.find (*it);
    if (ndit == dependencies_set.cend()) {
      dependencies_set.insert(*it);
      dependencies_no_dupl.insert(dependencies_no_dupl.begin(), *it);
    }
  }
  //print out the dependencies
  /*std::cout << "Formula: " << SBML_formulaToString(&formula) << std::endl ;
  for  (it=dependencies_no_dupl.begin(); it<dependencies_no_dupl.end(); it++) {
    std::cout << " " << *it ;
  }
  std::cout  <<  std::endl; */

  return dependencies_no_dupl;
}

} // end of namespace wcs

std::unordered_map<std::string, size_t> build_input_map(
  const std::string& formula,
  const std::set<std::string>& var_names)
{
  static const std::regex symbol_name("[\\w_]+"); // sequence of alnum or '_'
  std::unordered_map<std::string, size_t> input_map;
  std::sregex_token_iterator p(formula.cbegin(), formula.cend(), symbol_name);
  std::sregex_token_iterator e;
  size_t i = 0ul;
  for ( ; p!=e ; ++p) {
    for (const auto& vn : var_names) {
      if (vn == *p) {
        if (input_map.count(vn) == 0u) {
          input_map.insert(std::make_pair(vn, i++));
        }
        if (input_map.size() == var_names.size()) {
          return input_map;
        }
      }
    }
  }
  return input_map;
}

const std::string  wcs::generate_cxx_code::generate_code(const LIBSBML_CPP_NAMESPACE::Model& model)
{
  const ListOfCompartments* compartment_list = model.getListOfCompartments();
  const unsigned int compartmentsSize = compartment_list->size();
  const ListOfSpecies* species_list = model.getListOfSpecies();
  const unsigned int speciesSize = species_list->size();
  const ListOfReactions* reaction_list = model.getListOfReactions();
  const unsigned int reactionsSize = reaction_list->size();
  const ListOfParameters* parameter_list = model.getListOfParameters();
  const unsigned int parametersSize = parameter_list->size();
  const ListOfFunctionDefinitions* function_definition_list
    = model.getListOfFunctionDefinitions();
  const unsigned int functionsSize = function_definition_list->size();
  const ListOfRules* rules_list = model.getListOfRules();
  const unsigned int rulesSize = rules_list->size();
  const ListOfEvents* events_list = model.getListOfEvents();
  const unsigned int eventsSize = events_list->size();
  const ListOfInitialAssignments* assignments_list
    = model.getListOfInitialAssignments();
  const unsigned int assignmentsSize = assignments_list->size();

  const char * Real = generate_cxx_code::basetype_to_string<reaction_rate_t>::value;

  char * pointerFilename;
  pointerFilename = tmpnam (NULL);
  //std::cout << "name file:" << pointerFilename << ".cc" << "\n";
  std::stringstream ss;
  ss << pointerFilename;
  std::ofstream genfile;
  std::string filename = ss.str() + ".cc";
  genfile.open(filename);

  genfile << "/** Autogenerated source code, do not edit! \n */"
            << "\n\n//C++ includes\n"
            << "#include <vector>\n"
            << "#include <cmath>\n\n"
            << "#include <cstdio>\n\n"
            << "//Include the text of the SBML in here.\n"
            << "//Use \"xxd -i\" to get the text as a C symbol.\n"
            << "extern \"C\" \n{\n"
            << "  unsigned char __original_sbml[] = \"\";\n}\n\n"
            << "//Constants and other functions defined by the SBML standard.\n\n"
            << "//Get the correct floating point type from the code at runtime.\n"
            << "typedef " << Real << " reaction_rate_t;\n"
            //<< "typedef reaction_rate_t " << Real << ";\n\n"
            << "//Prototype all the functions\n";

  for (unsigned int ic = 0u; ic<functionsSize; ic++) {
    const LIBSBML_CPP_NAMESPACE::FunctionDefinition& function
      = *(function_definition_list->get(ic));
    genfile << "static inline " << Real << " " << function.getIdAttribute() <<"(";
    unsigned int num_arg = function.getNumArguments();

    for (unsigned int ia =0u; ia<num_arg; ia++) {
      genfile << Real << " " << function.getArgument(ia)->getName();
      if (ia!=num_arg-1) {
        genfile << ", ";
      }
    }
    genfile << ");\n";
  }

  using  model_constants = std::unordered_set<std::string>;
  typename model_constants::const_iterator constit;
  // A set for model constants
  model_constants sconstants;

  using  initial_assignments = std::unordered_map<std::string, const ASTNode *>;
  typename initial_assignments::const_iterator initassigit;
  // A map for initial_assignments
  initial_assignments sinitial_assignments;

  using  constant_init_ass = std::unordered_map<std::string, const ASTNode *>;
  typename constant_init_ass::const_iterator cinitasit;
  //  A map for constants in initial assignments
  constant_init_ass sconstant_init_assig;

  using  rate_rules = std::unordered_map <std::string, const ASTNode *>;
  typename rate_rules::const_iterator rrit;
  // A map for model rate rules
  rate_rules rate_rules_set;

  using  assignment_rules = wcs::map_symbol_to_ast_node_t;
  typename assignment_rules::const_iterator arit;
  // A map for model assignment rules
  assignment_rules assignment_rules_set;

  using  model_reactions = wcs::map_symbol_to_ast_node_t;
  typename model_reactions::const_iterator mrit;
  // A map for model reactions
  model_reactions model_reactions_set;

  std::unordered_set <std::string> good_params, used_params, model_species;
  std::unordered_set <std::string>::const_iterator upit, msit;


  // Put reactions in a map
  for (unsigned int ic = 0; ic < reactionsSize; ic++) {
    const LIBSBML_CPP_NAMESPACE::Reaction& reaction = *(reaction_list->get(ic));
    model_reactions_set.insert(std::make_pair(reaction.getIdAttribute(),
                                              reaction.getKineticLaw()->getMath()));

  }

  // Put initial assignments in map
  for (unsigned int ic = 0u; ic < assignmentsSize; ic++) {
    const LIBSBML_CPP_NAMESPACE::InitialAssignment& initialassignment
      = *(assignments_list->get(ic));
    sinitial_assignments.insert(std::make_pair(initialassignment.getSymbol(),
                                               initialassignment.getMath()));
  }

  // Put species in a set
  for (unsigned int ic = 0u; ic < speciesSize; ic++) {
    if (assignment_rules_set.find(species_list->get(ic)->getIdAttribute())
        == assignment_rules_set.cend())
    {
      model_species.insert(species_list->get(ic)->getIdAttribute());
    }
  }



  genfile << "\n";
  using  event_assignments = std::unordered_set<std::string>;
  typename event_assignments::const_iterator evassigit;
  // A set for model event assignments
  event_assignments m_ev_assig;

  //Put event assignments in a set
  for (unsigned int ic = 0u; ic < eventsSize; ic++) {
    const LIBSBML_CPP_NAMESPACE::Event& event = *(events_list->get(ic));
    const ListOfEventAssignments* eventassignments_list
      = event.getListOfEventAssignments();
    const unsigned int eventassignmentsSize = eventassignments_list->size();

    for (unsigned int ici = 0u; ici < eventassignmentsSize; ici++) {
      const LIBSBML_CPP_NAMESPACE::EventAssignment& eventassignment
        = *(eventassignments_list->get(ici));
      if (m_ev_assig.find(eventassignment.getVariable()) == m_ev_assig.cend()) {
        //genfile << eventassignment.getVariable() <<"\n";
        m_ev_assig.insert(eventassignment.getVariable());
      }
    }
  }


  // Put rate rules and assignements rules in maps
  for (unsigned int ic = 0; ic < rulesSize; ic++) {
    const LIBSBML_CPP_NAMESPACE::Rule& rule = *(rules_list->get(ic));
    if (rule.getType() == 0) { //rate_rule
      rate_rules_set.insert(std::make_pair(rule.getVariable(), rule.getMath()));
    }
    if (rule.getType() == 1) { //assignment_rule
      assignment_rules_set.insert(std::make_pair(rule.getVariable(), rule.getMath()));
    }

  }

  // Find used parameters in the rates
  for (unsigned int ic = 0u; ic < reactionsSize; ic++) {
    const LIBSBML_CPP_NAMESPACE::Reaction& reaction = *(reaction_list->get(ic));
    //genfile << "reaction" << reaction.getIdAttribute() <<": ";
    std::vector<std::string> dependencies_set
      = get_all_dependencies(*reaction.getKineticLaw()->getMath(),
                                  used_params,
                                  sinitial_assignments,
                                  assignment_rules_set,
                                  model_reactions_set);

    for (auto it = dependencies_set.crbegin();
         it != dependencies_set.crend(); ++it)
    {
      arit = assignment_rules_set.find(*it);
      upit = used_params.find(*it);
      rrit = rate_rules_set.find(*it);
      mrit = model_reactions_set.find(*it);
      msit = model_species.find(*it);
      evassigit = m_ev_assig.find(*it);
      if (arit == assignment_rules_set.cend() &&
          rrit == rate_rules_set.cend() &&
          mrit == model_reactions_set.cend() &&
          upit == used_params.cend() &&
          msit == model_species.cend() &&
          *it != "time" &&
          evassigit == m_ev_assig.cend())
      {
        used_params.insert(*it);
      }
    }
  }

  // Find used parameters in the differential rates
  for (unsigned int ic = 0; ic < rulesSize; ic++) {
    const LIBSBML_CPP_NAMESPACE::Rule& rule = *(rules_list->get(ic));
    if (rule.getType() == 0) { //rate_rule
      //genfile << "rule" << rule.getVariable() << " : " ;
      std::vector<std::string> dependencies_set
        = get_all_dependencies(*rule.getMath(),
                                    used_params,
                                    sinitial_assignments,
                                    assignment_rules_set,
                                    model_reactions_set);

      for (auto it = dependencies_set.crbegin();
           it!= dependencies_set.crend(); ++it)
      {
        arit = assignment_rules_set.find(*it);
        upit = used_params.find(*it);
        rrit = rate_rules_set.find(*it);
        mrit = model_reactions_set.find(*it);
        msit = model_species.find(*it);
        evassigit = m_ev_assig.find(*it);

        if (arit == assignment_rules_set.cend() &&
            rrit == rate_rules_set.cend() &&
            mrit == model_reactions_set.cend() &&
            upit == used_params.cend() &&
            msit == model_species.cend() &&
            *it != "time" &&
            evassigit == m_ev_assig.cend())
        {
          used_params.insert(*it);
        }
      }
    }
  }

  //genfile << "Size of used params: " << used_params.size() <<"\n" ;

  /*for (auto& x: used_params) {
    std::cout << x << "\n";
  }*/


  genfile << "//Define all the constants and initial states\n";
  for (unsigned int ic = 0u; ic < parametersSize; ic++) {
    const LIBSBML_CPP_NAMESPACE::Parameter& parameter
      = *(parameter_list->get(ic));

    if (parameter.getConstant()) {
      if (!isnan(parameter.getValue())) {  // Constants with value
        if (used_params.find(parameter.getIdAttribute()) !=
            used_params.cend())
        { //declare if it is used
          genfile << "static constexpr "<< Real << " " << parameter.getIdAttribute()
                    << " = " << parameter.getValue() << ";\n";
          good_params.insert(parameter.getIdAttribute());
        }
      } else {  //if no value for constants check initial assignments
        if (used_params.find(parameter.getIdAttribute()) !=
            used_params.cend())
        { //declare if it is used
          initassigit = sinitial_assignments.find(parameter.getIdAttribute());
          if (initassigit != sinitial_assignments.cend()) { //initial_assignments
            sconstant_init_assig.insert(std::make_pair(parameter.getIdAttribute(),
                                                       initassigit->second));
          }
        }
      }
    } else {
      if (!isnan(parameter.getValue())) {  // Variables with values
        genfile << "static constexpr " << Real << " _init_"
                  << parameter.getIdAttribute() << " = "
                  << parameter.getValue() << ";\n";
        good_params.insert(parameter.getIdAttribute());
      } else {  //if no value for these variables => they are assignment_rules
      }
    }
  }

  // Print out compartments initializations
  for (unsigned int ic = 0u; ic < compartmentsSize; ic++) {
    const LIBSBML_CPP_NAMESPACE::Compartment& compartment
      = *(compartment_list->get(ic));

    if (rate_rules_set.find(compartment.getIdAttribute()) ==
        rate_rules_set.cend())
    {
      genfile << "static constexpr " << Real << " " << compartment.getIdAttribute()
                << " = " << compartment.getSize() << ";\n";
      good_params.insert(compartment.getIdAttribute());
    } else { //if it is defined in a rate_rule
      genfile << "static constexpr " << Real << " _init_"
                << compartment.getIdAttribute() << " = "
                << compartment.getSize() << ";\n";
      good_params.insert(compartment.getIdAttribute());
    }
  }


  /// Print out Species initializations
  for (unsigned int ic = 0u; ic < speciesSize; ic++) {
    const LIBSBML_CPP_NAMESPACE::Species& species = *(species_list->get(ic));

    if (species.getInitialAmount() == 0) {
      genfile << "static constexpr " << Real  << " _init_"
                << species.getIdAttribute() << " = 0;\n";
      //good_params.insert(species.getIdAttribute());
    } else {
      if (!isnan(species.getInitialAmount())) {
        // First way for printing the value
        std::stringstream ss;
        ss << species.getInitialAmount();
        std::string speciesvalue = ss.str();

        // Second way. Print the value with big precision
        //char svalue [100];
        //char * spname = new char [species.getIdAttribute().length()+1];
        //strcpy (spname, species.getIdAttribute().c_str());
        //char * compname = new char [species.getCompartment().length()+1];
        //strcpy (compname, species.getCompartment().c_str());

        //sprintf (svalue, "%.16f", species.getInitialAmount());
        //printf ("static constexpr Real _init_%s = %s/%s; ", spname, svalue, compname);*/

        // Third way
        //std::ostringstream out;
        //out.precision(15);
        //out << std::fixed << species.getInitialAmount();

        genfile << "static constexpr " << Real << " _init_"
                  << species.getIdAttribute() << " = "
                  << speciesvalue << "/"
                  //<< svalue << "/"
                  //<< out.str() << "/"
                  << species.getCompartment() << ";\n";



      } else if (!isnan(species.getInitialConcentration())) {
        if (species.getInitialConcentration() == 0) {
          genfile << "static constexpr " << Real << " _init_"
                    << species.getIdAttribute() << " = 0;\n";
        } else {
          // First way for printing the value
          std::stringstream ss;
          ss << species.getInitialConcentration();
          std::string speciesvalue = ss.str();

          genfile << "static constexpr " << Real << " _init_"
                    << species.getIdAttribute() << " = "
                    << speciesvalue << "/" << species.getCompartment() << ";\n";
        }
      }
    }
  }

  /** for (auto x: sconstant_init_assig) {
      genfile << x.first << "\n";
  }*/


  // A set with parameters that they have to be deleted after the for loop
  std::unordered_set <std::string> param_to_del;

  // Print out constants in initial assignemnts
  for (auto x = sconstant_init_assig.begin(); x != sconstant_init_assig.end(); ++x){
    if (param_to_del.find(x->first) == param_to_del.cend()) {
      if (rate_rules_set.find(SBML_formulaToString(x->second)) ==
          rate_rules_set.cend())
      {
        //if (used_params.find(SBML_formulaToString(x->second))
        //!= used_params.cend()){ //declare if it is used
        std::vector<std::string> dependencies_set
          = get_all_dependencies(*x->second,
                                      good_params,
                                      sconstant_init_assig,
                                      assignment_rules_set,
                                      model_reactions_set);

        for (auto it = dependencies_set.crbegin();
             it != dependencies_set.crend(); ++it)
        {
          arit = assignment_rules_set.find(*it);
          cinitasit = sconstant_init_assig.find(*it);
          if (arit != assignment_rules_set.cend()) {
            /**const std::string toReplace("species");
               // for the replacement of "species" label
            size_t pos = formula.find(toReplace);
            while (pos < formula.length()) {
              if (pos != std::string::npos) {
                formula.replace(pos,toReplace.length(),"");
                pos = formula.find(toReplace);
              }
            }*/
            genfile << "static constexpr " << Real << " " << arit->first << " = "
                      << SBML_formulaToString(arit->second) << ";\n";
          }
          if (cinitasit != sconstant_init_assig.cend()){
            genfile << "static constexpr " << Real << " "
                      << cinitasit->first << " = "
                      << SBML_formulaToString(cinitasit->second) << ";\n";
            good_params.insert(cinitasit->first);
            param_to_del.insert(cinitasit->first);
          }
        }
        genfile << "static constexpr " << Real << " " << x->first << " = "
                  << SBML_formulaToString(x->second) << ";\n";
        good_params.insert(x->first);
        param_to_del.insert(x->first);
      //}
      } else {
        genfile << "static constexpr " << Real << " " << x->first << " = _init_"
                  << SBML_formulaToString(x->second) << ";\n";
        good_params.insert(x->first);
        param_to_del.insert(x->first);
      }
    }
  }

  // Remove parameters that have been declared
  for (auto& x: param_to_del) {
    sconstant_init_assig.erase(x);
  }

  genfile << "\n";
  genfile << "extern \"C\" \n{\n";
  // Print out  alias for event assignments
  for  (const std::string& x: m_ev_assig) {
    std::string e_label = x;
    genfile << "  " << Real << " __init_" << e_label << " = _init_"
              << e_label <<";\n";
    // we delete these from good_params because we have to declare in rates
    good_params.erase(e_label);
  }

  // Print out alias for rate rules
  for  (auto& x: rate_rules_set) {
    std::string r_label = x.first;
    genfile << "  " << Real << " __init_" << r_label << " = _init_"
              << r_label <<";\n";
    // we delete these from good_params because we have to declare in rates
    good_params.erase(r_label);
  }

  // Print out alias for species thay are not defined in assignment rules
  for (unsigned int ic = 0u; ic < speciesSize; ic++) {
    if (assignment_rules_set.find(species_list->get(ic)->getIdAttribute()) ==
        assignment_rules_set.cend())
    {
        genfile << "  " << Real << " __init_"
                  << species_list->get(ic)->getIdAttribute() << " = _init_"
                  << species_list->get(ic)->getIdAttribute() << ";\n";
    }
  }
  genfile << "}";

  genfile << "\n//Define the functions\n";
  for (unsigned int ic = 0u; ic < functionsSize; ic++) {
    const LIBSBML_CPP_NAMESPACE::FunctionDefinition& function
      = *(function_definition_list->get(ic));

    genfile << "static inline " << Real << " " << function.getIdAttribute() << "(";
    unsigned int num_arg = function.getNumArguments();

    for (unsigned int ia = 0u; ia < num_arg; ia++) {
      genfile << Real << " " << function.getArgument(ia)->getName();
      if (ia != num_arg-1) {
        genfile << ", ";
      }
    }

    genfile << ") {\n";
    genfile << "  return " << SBML_formulaToString(function.getBody()) << ";";
    genfile << "\n}\n\n";
  }


  /**genfile << "Size of used params: " << used_params.size() <<"\n" ;

  for (auto& x: good_params) {
    if (used_params.find(x) == used_params.cend()){
      std::cout << x << "\n";
    }
  }
  genfile << "Size of good params: " << good_params.size() <<"\n\n" ;*/


  genfile << "//Define the rates\n";
  for (unsigned int ic = 0u; ic < reactionsSize; ic++) {
    const LIBSBML_CPP_NAMESPACE::Reaction& reaction = *(reaction_list->get(ic));

    genfile << "extern \"C\" " << Real << " wcs__rate_" << reaction.getIdAttribute()
              << "(const std::vector<" << Real <<">& __input) {\n";
    genfile << "  int __ii=0;\n";
    //genfile << "printf(\"nfunction: %s \\n\", \"" << reaction.getIdAttribute() << "\");\n";
    std::vector<std::string> dependencies_set
      = get_all_dependencies(*reaction.getKineticLaw()->getMath(),
                                  good_params,
                                  sconstant_init_assig,
                                  assignment_rules_set,
                                  model_reactions_set);

    std::set<std::string> var_names;
    // print the function parameters with the order met in the formula
    //put first to a set all variables taking an input
    for (auto it = dependencies_set.crbegin();
         it != dependencies_set.crend(); ++it)
    {
      arit = assignment_rules_set.find(*it);
      if (arit == assignment_rules_set.cend()){
        var_names.insert(*it);
      }
    }
    std::string formula = SBML_formulaToString(reaction.getKineticLaw()->getMath());
    const auto input_map = build_input_map(formula, var_names);

    std::vector<std::string> var_names_ord (input_map.size());
    for (auto& x: input_map) {
      var_names_ord[x.second] = x.first;
    }

    //print first parameters in formula taking input
    std::vector<std::string>::iterator it;
    for (it = var_names_ord.begin(); it < var_names_ord.end(); it++){
      genfile << "  " << Real << " " << *it <<" = __input[__ii++];\n";
      var_names.erase(*it);
    }
    //print rest parameters (not in formula) taking input
    for (auto& x: var_names) {
      genfile << "  " << Real << " " << x <<" = __input[__ii++];\n";
    }

    //print the rest of the parameters not taking an input
    for (auto it = dependencies_set.crbegin();
         it != dependencies_set.crend(); ++it)
    {
      arit = assignment_rules_set.find(*it);
      if (arit != assignment_rules_set.cend()){
        genfile << "  " << Real << " " << arit->first << " = "
                  << SBML_formulaToString(arit->second) << ";\n";
      }
    }

    genfile << "  " << Real << " " << reaction.getIdAttribute() << " = "
              << SBML_formulaToString(reaction.getKineticLaw()->getMath())
              << ";\n";
    genfile << "  return " << reaction.getIdAttribute() << ";\n" << "}\n\n";
  }

  genfile << "\n//Define differential rates\n";
  for (unsigned int ic = 0u; ic < rulesSize; ic++) {
    const LIBSBML_CPP_NAMESPACE::Rule& rule = *(rules_list->get(ic));
    if (rule.getType()==0) { //rate_rule
      genfile << "extern \"C\" " << Real << " wcs__rate_" << rule.getVariable()
                << "(const std::vector<" << Real << ">& __input) {\n";
      genfile << "  int __ii=0;\n";

      std::vector<std::string> dependencies_set
        = get_all_dependencies(*rule.getMath(),
                                    good_params,
                                    sconstant_init_assig,
                                    assignment_rules_set,
                                    model_reactions_set);
      std::set<std::string> var_names;
      // print the function parameters with the order met in the formula
      // put first to a set all variables taking an input
      for (auto it = dependencies_set.crbegin();
           it!= dependencies_set.crend(); ++it)
      {
        arit = assignment_rules_set.find(*it);
        mrit = model_reactions_set.find(*it);
        if (arit == assignment_rules_set.cend() && mrit == model_reactions_set.cend()){
          var_names.insert(*it);
        }
      }
      std::string formula = SBML_formulaToString(rule.getMath());
      const auto input_map = build_input_map(formula, var_names);

      std::vector<std::string> var_names_ord (input_map.size());
      for (auto& x: input_map) {
        var_names_ord[x.second] = x.first;
      }

      //print first parameters in formula taking input
      std::vector<std::string>::iterator it;
      for (it = var_names_ord.begin(); it < var_names_ord.end(); it++){
        genfile << "  " << Real << " " << *it <<" = __input[__ii++];\n";
        var_names.erase(*it);
      }
      //print rest parameters (not in formula) taking input
      for (auto& x: var_names) {
        genfile << "  " << Real << " " << x <<" = __input[__ii++];\n";
      }


      // Print the rest of the parameters not taking an input
      for (auto it = dependencies_set.crbegin();
           it!= dependencies_set.crend(); ++it)
      {
        arit = assignment_rules_set.find(*it);
        mrit = model_reactions_set.find(*it);
        if (arit != assignment_rules_set.cend()){
          genfile << "  " << Real << " " << arit->first << " = "
                    << SBML_formulaToString(arit->second) << ";\n";
        } else if (mrit != model_reactions_set.cend()){
          genfile << "  " << Real << " " << mrit->first << " = "
                    << SBML_formulaToString(mrit->second) << ";\n";
        }
        /*else {
          genfile << "  " << Real << " " << *it <<" = __input[__ii++];\n";
        }*/
      }

      genfile << "\n  return " << SBML_formulaToString(rule.getMath())
                << ";\n" << "}\n\n";
    }
  }

  genfile.close();
  return filename;

}

const std::string  wcs::generate_cxx_code::compile_code(const std::string generated_filename)
{
  std::string dir;
  std::string stem;
  std::string ext;
  extract_file_component(generated_filename, dir, stem, ext);

  //std::string lib_filename1 = dir + stem + ".o";
  //std::string lib_filename2 = dir + stem + ".so";

  std::string lib_filename1 = stem + ".o";
  std::string lib_filename2 = stem + ".so";

  std::string system1 = std::string(CMAKE_CXX_COMPILER) + " -std=c++17 -g -fPIC " + CMAKE_CXX_SHARED_LIBRARY_FLAGS + " -c " + generated_filename;
  std::string system2 = std::string(CMAKE_CXX_COMPILER) + " -std=c++17 -g -fPIC " + CMAKE_CXX_SHARED_LIBRARY_FLAGS  + " -shared -export-dynamic " + lib_filename1 + " -o " + lib_filename2;

  system(system1.c_str());
  system(system2.c_str());

  return lib_filename2;
}



#endif // defined(WCS_HAS_SBML)
