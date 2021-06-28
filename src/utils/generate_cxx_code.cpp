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
#include <algorithm> // std::min
#include <cstdlib> // system, realpath
#include <climits> // PATH_MAX
#include <cstring> // strncpy
#include <fcntl.h> // O_CREAT
#include <unistd.h> // close
#include <sys/wait.h> // WEXITSTATUS
#include <set>
#include <regex>
#include <string>
#include "wcs_types.hpp"
#include "sbml_functions.hpp"
#include <thread> // std::thread::hardware_concurrency

#if defined(WCS_HAS_SBML)

LIBSBML_CPP_NAMESPACE_USE

namespace wcs {

using map_symbol_to_ast_node_t
  = std::unordered_map <std::string, const LIBSBML_CPP_NAMESPACE::ASTNode *>;
using params_map_t = std::unordered_map <std::string, std::vector<std::string>>;
using rate_rules_dep_t = std::unordered_map <std::string, std::set<std::string>>;
using event_assignments_t = std::unordered_set<std::string>;
using initial_assignments_t = wcs::map_symbol_to_ast_node_t;
using assignment_rules_t = wcs::map_symbol_to_ast_node_t;
using model_reactions_t = wcs::map_symbol_to_ast_node_t;
using rate_rules_t = wcs::map_symbol_to_ast_node_t;
using constant_init_ass_t = map_symbol_to_ast_node_t;

template<>
const char* generate_cxx_code::basetype_to_string<double>::value = "double";
template<>
const char* generate_cxx_code::basetype_to_string<float>::value = "float";
typedef reaction_rate_t ( * rate_function_pointer)(const std::vector<reaction_rate_t>&);

void
generate_cxx_code::get_dependencies(
  const LIBSBML_CPP_NAMESPACE::ASTNode& math,
  std::vector<std::string> & math_elements,
  const std::unordered_set<std::string> & good,
  const map_symbol_to_ast_node_t & constant_init_assig,
  const map_symbol_to_ast_node_t & variables_assig_rul,
  const map_symbol_to_ast_node_t & model_reactions_s,
  const std::unordered_set<std::string> & local_params)
{
  if (math.getNumChildren() >= 1u) {
    for (unsigned int n = 0u; n < math.getNumChildren(); ++n) {
      const LIBSBML_CPP_NAMESPACE::ASTNode*  math_c = math.getChild(n);
      get_dependencies(*math_c,
                       math_elements,
                       good,
                       constant_init_assig,
                       variables_assig_rul,
                       model_reactions_s,
                       local_params);
    }
  } else {
    if (!math.isNumber()) {
      const char* formula = SBML_formulaToString(&math);
      map_symbol_to_ast_node_t::const_iterator inasit, asruit, mreactit;
      inasit = constant_init_assig.find(formula);
      asruit = variables_assig_rul.find(formula);
      mreactit = model_reactions_s.find(formula);
      std::unordered_set<std::string>::const_iterator git = good.find(formula);
      std::unordered_set<std::string>::const_iterator lpit = local_params.find(formula);
      //std::unordered_set<std::string> no_duplicates_vector;
      if (git == good.cend()) {
        if (inasit != constant_init_assig.cend()) {
          math_elements.insert(math_elements.cbegin(), formula);
          get_dependencies(*inasit->second,
                           math_elements,
                           good,
                           constant_init_assig,
                           variables_assig_rul,
                           model_reactions_s,
                           local_params);

        } else if (asruit != variables_assig_rul.cend()) {
          math_elements.insert(math_elements.cbegin(), formula);
          get_dependencies(*asruit->second,
                           math_elements,
                           good,
                           constant_init_assig,
                           variables_assig_rul,
                           model_reactions_s,
                           local_params);

        } else if (mreactit != model_reactions_s.cend()) {
          math_elements.insert(math_elements.cbegin(), formula);
          get_dependencies(*mreactit->second,
                           math_elements,
                           good,
                           constant_init_assig,
                           variables_assig_rul,
                           model_reactions_s,
                           local_params);

        } else if (lpit != local_params.cend()) {

        } else {
          math_elements.insert(math_elements.cbegin(), formula);
        }
      }
    }
  }
}


std::vector<std::string>
generate_cxx_code::get_all_dependencies(
  const LIBSBML_CPP_NAMESPACE::ASTNode& formula,
  const std::unordered_set<std::string> & good,
  const map_symbol_to_ast_node_t & constant_init_assig,
  const map_symbol_to_ast_node_t & variables_assig_rul,
  const map_symbol_to_ast_node_t & model_reactions_s,
  const std::unordered_set<std::string> & local_params)
{

  std::vector<std::string> dependencies, dependencies_no_dupl;
  std::unordered_set<std::string> dependencies_set;
  std::unordered_set<std::string>::const_iterator ndit;

  get_dependencies(formula,
                   dependencies,
                   good,
                   constant_init_assig,
                   variables_assig_rul,
                   model_reactions_s,
                   local_params);

  std::vector<std::string>::const_iterator it;
  //remove duplicates
  for  (it = dependencies.cbegin(); it < dependencies.cend(); it++) {
    ndit = dependencies_set.find (*it);
    if (ndit == dependencies_set.cend()) {
      dependencies_set.insert(*it);
      dependencies_no_dupl.insert(dependencies_no_dupl.cbegin(), *it);
    }
  }
  //print out the dependencies
  /*std::cout << "Formula: " << SBML_formulaToString(&formula) << std::endl ;
  for  (it=dependencies_no_dupl.cbegin(); it<dependencies_no_dupl.cend(); it++) {
    std::cout << " " << *it ;
  }
  std::cout  <<  std::endl; */

  return dependencies_no_dupl;
}

static double count2Concentration_converter (
  const LIBSBML_CPP_NAMESPACE::Model& model,
  const LIBSBML_CPP_NAMESPACE::Species& species,
  std::string& compartment,
  double& comp_unit)
{
  const LIBSBML_CPP_NAMESPACE::ListOfCompartments* compartment_list
   = model.getListOfCompartments();
  const LIBSBML_CPP_NAMESPACE::ListOfUnitDefinitions* unit_definition_list
   = model.getListOfUnitDefinitions(); 
  compartment = species.getCompartment();
  double compartment_size = compartment_list->get(species.getCompartment())->getSize();
  double avog_num = 6.02214e+23;
  
               
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
  // double comp_unit = 1.0;
  comp_unit = 1.0;
  for (unsigned int iu = 0u; iu < unitsSize; iu++) { 
    const LIBSBML_CPP_NAMESPACE::Unit* unit = unit_list->get(iu);
    comp_unit = comp_unit * pow(unit->getMultiplier()*pow(10,unit->getScale()),unit->getExponent());
  }  
  return avog_num * compartment_size * comp_unit; 
}


void
generate_cxx_code::return_denominators(
  const LIBSBML_CPP_NAMESPACE::ASTNode& formula,
  std::vector<std::string> & denominators,
  const std::string& reaction_name)
{

  if (formula.getNumChildren() >= 1u) {
    if (formula.getType() == AST_DIVIDE) {
      //std::cout << " " << SBML_formulaToString(formula.getChild(1)) << std::endl;
      const LIBSBML_CPP_NAMESPACE::ASTNode*  denominator = formula.getChild(1);
      if (denominator->isNumber()) {
        if (denominator->getValue() == 0) {
          throw std::overflow_error("Divide by zero exception in reaction " + reaction_name);
        }
      } else {
        if ( std::find(denominators.cbegin(), denominators.cend(),
        SBML_formulaToString(denominator)) == denominators.cend()) {
          std::string new_denominator = sbml_functions::SBML_formulaToString_wcs(denominator);
          denominators.push_back(new_denominator);
        }
      }
    }
    for (unsigned int n = 0u; n < formula.getNumChildren(); ++n) {
      const LIBSBML_CPP_NAMESPACE::ASTNode*  math_c = formula.getChild(n);
      return_denominators(*math_c, denominators, reaction_name);
    }
  }
}

std::vector<std::string>
generate_cxx_code::return_all_denominators(
  const LIBSBML_CPP_NAMESPACE::ASTNode& formula,
  const std::string& reaction_name)
{
  std::vector<std::string> denominators;
  return_denominators(formula, denominators, reaction_name);
  return denominators;
}

static std::unordered_map<std::string, size_t> build_input_map(
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

struct strcomp {
  bool operator() (const std::string& lhs, const std::string& rhs) const
  {return lhs.size()>rhs.size();}
};

void read_ast_node(
  const LIBSBML_CPP_NAMESPACE::ASTNode& math,
  std::set<std::string, strcomp>& math_elements)
{
  if (math.getNumChildren() >= 1u) {
    for (unsigned int n = 0u; n < math.getNumChildren(); ++n) {
      const LIBSBML_CPP_NAMESPACE::ASTNode* math_c = math.getChild(n);
      read_ast_node(*math_c, math_elements);
    }
  } else {
    if (!math.isNumber()) {
      char* formula = SBML_formulaToString(&math);
      std::set<std::string,strcomp>::const_iterator mit
        = math_elements.find (formula);
      if (mit == math_elements.cend()) {
        math_elements.insert(formula);
      }
    }
  }
}

// Update with the correct scope a math formula
static void update_scope_ast_node(
  LIBSBML_CPP_NAMESPACE::ASTNode& math,
  const std::unordered_set<std::string>& wcs_var,
  const std::unordered_set<std::string>& wcs_const,
  const std::unordered_set<std::string>& local_params)
{
  if (math.getNumChildren() >= 1u) {
    for (unsigned int n = 0u; n < math.getNumChildren(); ++n) {
      LIBSBML_CPP_NAMESPACE::ASTNode* math_c = math.getChild(n);
      update_scope_ast_node(*math_c, wcs_var, wcs_const, local_params);
    }
  } else {
    if (!math.isNumber()) {
      char* formula = SBML_formulaToString(&math);
      std::unordered_set<std::string>::const_iterator wcs_var_it, wcs_const_it, local_params_it;
      wcs_var_it = wcs_var.find(formula);
      wcs_const_it = wcs_const.find(formula);
      local_params_it = local_params.find(formula);
      const char* nodename = math.getName();
      if (local_params_it == local_params.cend()) {
        if (wcs_const_it != wcs_const.cend()) {
          std::string new_nodename = "WCS_GLOBAL_CONST::";
          new_nodename += nodename;
          math.setName(new_nodename.c_str());
        }
        if (wcs_var_it != wcs_var.cend()) {

          std::string new_nodename = "wcs_global_var.";
          new_nodename += nodename;
          math.setName(new_nodename.c_str());
        }
      }
    }

  }
}

// Update with the correct scope a element string
static void update_scope_str(
  std::string& math,
  const std::unordered_set<std::string>& wcs_var,
  const std::unordered_set<std::string>& wcs_const,
  const std::unordered_set<std::string>& local_params)
{
  std::unordered_set<std::string>::const_iterator wcs_var_it, wcs_const_it, local_params_it;
  wcs_var_it = wcs_var.find(math);
  wcs_const_it = wcs_const.find(math);
  local_params_it = local_params.find(math);
  if (local_params_it == local_params.cend()) {
    if (wcs_const_it != wcs_const.cend()) {
      math = "WCS_GLOBAL_CONST::" + math;
    }
    if (wcs_var_it != wcs_var.cend()) {
      math = "wcs_global_var." + math;
    }
  }
}

// Include init for initializations with rate rules
static void include_init_for_rate_rules(
  LIBSBML_CPP_NAMESPACE::ASTNode& math,
  const std::unordered_map <std::string, const ASTNode *>& raterules)
{
  if (math.getNumChildren() >= 1u) {
    for (unsigned int n = 0u; n < math.getNumChildren(); ++n) {
      LIBSBML_CPP_NAMESPACE::ASTNode* math_c = math.getChild(n);
      include_init_for_rate_rules(*math_c, raterules);
    }
  } else {
    if (!math.isNumber()) {
      char* formula = SBML_formulaToString(&math);
      std::unordered_map <std::string, const ASTNode *>::const_iterator rrit;
      rrit = raterules.find(formula);
      const char* nodename = math.getName();
      if (rrit != raterules.cend()) {
        std::string new_nodename = "_init_";
        new_nodename += nodename;
        math.setName(new_nodename.c_str());
      }
    }
  }
}

void generate_cxx_code::find_used_params(
  const LIBSBML_CPP_NAMESPACE::Model& model,
  std::unordered_set <std::string>& used_params,
  const initial_assignments_t& sinitial_assignments,
  const assignment_rules_t& assignment_rules_map,
  const model_reactions_t& model_reactions_map,
  const rate_rules_t& rate_rules_map,
  const std::unordered_set <std::string>& model_species,
  const event_assignments_t& ev_assign)
{
  const ListOfReactions* reaction_list = model.getListOfReactions();
  const unsigned int num_reactions = reaction_list->size();
  const ListOfRules* rules_list = model.getListOfRules();
  const unsigned int num_rules = rules_list->size();
  typename assignment_rules_t::const_iterator arit;
  typename model_reactions_t::const_iterator mrit;
  std::unordered_set <std::string>::const_iterator upit, msit;
  typename rate_rules_t::const_iterator rrit;
  typename event_assignments_t::const_iterator evassigit;
  // Find used parameters in the rates
  for (unsigned int ic = 0u; ic < num_reactions; ic++) {
    const LIBSBML_CPP_NAMESPACE::Reaction& reaction = *(reaction_list->get(ic));
    // const LIBSBML_CPP_NAMESPACE::ListOfLocalParameters* local_parameter_list
    //   = reaction.getKineticLaw()->getListOfLocalParameters();
    // unsigned int num_localparameters = local_parameter_list->size();
    // //genfile << "reaction" << reaction.getIdAttribute() <<": ";
    using reaction_parameters_t = std::unordered_set<std::string>;
    reaction_parameters_t lpset;
    // for (unsigned int pi = 0u; pi < num_localparameters; pi++) {
    //   const LIBSBML_CPP_NAMESPACE::LocalParameter* localparameter = local_parameter_list->get(pi);
    //   lpset.insert(localparameter->getIdAttribute());
    // Check SBML level for local parameters
    if (model.getLevel() > 2) {
      const LIBSBML_CPP_NAMESPACE::ListOfLocalParameters* local_parameter_list
        = reaction.getKineticLaw()->getListOfLocalParameters();
      unsigned int num_localparameters = local_parameter_list->size();
      //genfile << "reaction" << reaction.getIdAttribute() <<": ";
      
      for (unsigned int pi = 0u; pi < num_localparameters; pi++) {
        const LIBSBML_CPP_NAMESPACE::LocalParameter* localparameter = local_parameter_list->get(pi);
        lpset.insert(localparameter->getIdAttribute());
      }
    } else {
      const LIBSBML_CPP_NAMESPACE::ListOfParameters* local_parameter_list
        = reaction.getKineticLaw()->getListOfParameters();
      unsigned int num_localparameters = local_parameter_list->size();
      //genfile << "reaction" << reaction.getIdAttribute() <<": ";
      
      for (unsigned int pi = 0u; pi < num_localparameters; pi++) {
        const LIBSBML_CPP_NAMESPACE::Parameter* localparameter = local_parameter_list->get(pi);
        lpset.insert(localparameter->getIdAttribute());
      } 
    }
    std::vector<std::string> dependencies_set
      = get_all_dependencies(*reaction.getKineticLaw()->getMath(),
                                  used_params,
                                  sinitial_assignments,
                                  assignment_rules_map,
                                  model_reactions_map,
                                  lpset);

    for (auto it = dependencies_set.crbegin();
         it != dependencies_set.crend(); ++it)
    {
      arit = assignment_rules_map.find(*it);
      upit = used_params.find(*it);
      rrit = rate_rules_map.find(*it);
      mrit = model_reactions_map.find(*it);
      msit = model_species.find(*it);
      evassigit = ev_assign.find(*it);
      if (arit == assignment_rules_map.cend() &&
          rrit == rate_rules_map.cend() &&
          mrit == model_reactions_map.cend() &&
          upit == used_params.cend() &&
          msit == model_species.cend() &&
          //*it != "time" &&   //uncomment for events
          evassigit == ev_assign.cend())
      {
        used_params.insert(*it);
      }
    }
  }

  // Find used parameters in the differential rates
  for (unsigned int ic = 0; ic < num_rules; ic++) {
    const LIBSBML_CPP_NAMESPACE::Rule& rule = *(rules_list->get(ic));
    if (rule.getType() == 0) { //rate_rule
      //genfile << "rule" << rule.getVariable() << " : " ;
      std::vector<std::string> dependencies_set
        = get_all_dependencies(*rule.getMath(),
                                    used_params,
                                    sinitial_assignments,
                                    assignment_rules_map,
                                    model_reactions_map,
                                    {});

      for (auto it = dependencies_set.crbegin();
           it!= dependencies_set.crend(); ++it)
      {
        arit = assignment_rules_map.find(*it);
        upit = used_params.find(*it);
        rrit = rate_rules_map.find(*it);
        mrit = model_reactions_map.find(*it);
        msit = model_species.find(*it);
        evassigit = ev_assign.find(*it);

        if (arit == assignment_rules_map.cend() &&
            rrit == rate_rules_map.cend() &&
            mrit == model_reactions_map.cend() &&
            upit == used_params.cend() &&
            msit == model_species.cend() &&
            //*it != "time" &&  //uncomment for events
            evassigit == ev_assign.cend())
        {
          used_params.insert(*it);
        }
      }
    }
  }

}

void generate_cxx_code::print_constants_and_initial_states(
  const LIBSBML_CPP_NAMESPACE::Model& model,
  std::ostream & genfile_hdr,
  std::ostream & genfile_impl,
  map_symbol_to_ast_node_t & sconstant_init_assig,
  const initial_assignments_t& sinitial_assignments,
  const assignment_rules_t& assignment_rules_map,
  const std::unordered_set <std::string>& used_params,
  std::unordered_set <std::string>& good_params,
  const model_reactions_t & model_reactions_map,
  const rate_rules_t& rate_rules_map,
  std::unordered_set<std::string>& wcs_all_const,
  std::unordered_set<std::string>& wcs_all_var,
  const event_assignments_t& ev_assign)
{
  const char* Real = generate_cxx_code::basetype_to_string<reaction_rate_t>::value;
  const ListOfParameters* parameter_list = model.getListOfParameters();
  const unsigned int num_parameters = parameter_list->size();
  const ListOfCompartments* compartment_list = model.getListOfCompartments();
  const unsigned int num_compartments = compartment_list->size();
  const ListOfEvents* events_list = model.getListOfEvents();
  const unsigned int num_events = events_list->size();
  typename assignment_rules_t::const_iterator arit;
  typename event_assignments_t::const_iterator evassigit;
  typename initial_assignments_t::const_iterator initassigit;
  typename constant_init_ass_t::const_iterator cinitasit;
  std::unordered_set <std::string> no_used_params;
  constant_init_ass_t sconstant_init_assig_notused;


  //Vector for const and var to keep the right order
  std::vector<std::string> wcs_const_exp, wcs_var_exp;

  // Maps for model constants with initial value and model variables with initial value
  using  wcs_struct_map_val_t = std::unordered_map <std::string,double>;
  typename wcs_struct_map_val_t::const_iterator wcs_const_itv, wcs_var_itv;
  std::unordered_map <std::string, double> wcs_const, wcs_var;

  // Maps for model constants and model variables with initial value other parameters
  using  wcs_struct_map_t = std::unordered_map <std::string, const ASTNode *>;
  typename wcs_struct_map_t::const_iterator wcs_const_it, wcs_var_it;
  wcs_struct_map_t wcs_const_exp_map, wcs_var_exp_map;

  for (unsigned int ic = 0u; ic < num_parameters; ic++) {
    const LIBSBML_CPP_NAMESPACE::Parameter& parameter
      = *(parameter_list->get(ic));

    if (parameter.getConstant()) {
      if (!isnan(parameter.getValue())) {  // Constants with value
        if (used_params.find(parameter.getIdAttribute()) !=
            used_params.cend())
        { //declare if it is used
          //std::ostringstream streamObj2;
          // Set Fixed -Point Notation
          //streamObj2 << std::fixed;
          //Add double to stream
          //streamObj2 << parameter.getValue();

          wcs_const.insert(std::make_pair(parameter.getIdAttribute(), parameter.getValue()));
          good_params.insert(parameter.getIdAttribute());
        } else {
          no_used_params.insert(parameter.getIdAttribute());
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
        } else {
          initassigit = sinitial_assignments.find(parameter.getIdAttribute());
          sconstant_init_assig_notused.insert(std::make_pair(parameter.getIdAttribute(),
                                                       initassigit->second));
          no_used_params.insert(parameter.getIdAttribute());
        }
      }
    } else {
      if (!isnan(parameter.getValue())) {  // Variables with values
        wcs_const.insert(std::make_pair("_init_" + parameter.getIdAttribute(),
        parameter.getValue()));
      } else {  //if no value for these variables => they are assignment_rules
      }
    }
  }

  // Compartments initializations
  for (unsigned int ic = 0u; ic < num_compartments; ic++) {
    const LIBSBML_CPP_NAMESPACE::Compartment& compartment
      = *(compartment_list->get(ic));

    if (rate_rules_map.find(compartment.getIdAttribute()) ==
        rate_rules_map.cend())
    {
      good_params.insert(compartment.getIdAttribute());
      wcs_const.insert(std::make_pair(compartment.getIdAttribute(), compartment.getSize()));
    } else { //if it is defined in a rate_rule
      wcs_const.insert(std::make_pair(" _init_" + compartment.getIdAttribute(),
      compartment.getSize()));
    }
  }

  // A set with parameters that they have to be deleted after the for loop
  std::unordered_set <std::string> param_to_del;

  // Constants in initial assignemnts
  for (auto x = sconstant_init_assig.cbegin(); x != sconstant_init_assig.cend(); ++x){
    if (param_to_del.find(x->first) == param_to_del.cend()) {
      if (rate_rules_map.find(SBML_formulaToString(x->second)) ==
          rate_rules_map.cend())
      {
        std::vector<std::string> dependencies_set
          = get_all_dependencies(*x->second,
                                      good_params,
                                      sconstant_init_assig,
                                      assignment_rules_map,
                                      model_reactions_map,
                                      {});

        for (auto it = dependencies_set.crbegin();
             it != dependencies_set.crend(); ++it)
        {
          arit = assignment_rules_map.find(*it);
          cinitasit = sconstant_init_assig.find(*it);
          if (arit != assignment_rules_map.cend()) {
            wcs_var_exp.push_back(arit->first);
            wcs_var_exp_map.insert(std::make_pair(arit->first,arit->second));
          }
          if (cinitasit != sconstant_init_assig.cend()){
            wcs_const_exp.push_back(cinitasit->first);
            wcs_const_exp_map.insert(std::make_pair(cinitasit->first,cinitasit->second));
            good_params.insert(cinitasit->first);
            param_to_del.insert(cinitasit->first);
          }
        }
        wcs_const_exp.push_back(x->first);
        wcs_const_exp_map.insert(std::make_pair(x->first,x->second));
        good_params.insert(x->first);
        param_to_del.insert(x->first);
      } else {
        wcs_const_exp.push_back(x->first);
        wcs_const_exp_map.insert(std::make_pair(x->first, x->second));
        good_params.insert(x->first);
        param_to_del.insert(x->first);
      }
    }
  }

  // No used parameters
  /*for (auto& x: no_used_params) {
    genfile << "\n" << x <<"," ;
  }*/

  // print no used parameters which are only used in the events formula
  if ( ev_assign.size() > 0) {
    for (unsigned int ic = 0u; ic < num_events; ic++) {
      const LIBSBML_CPP_NAMESPACE::Event& event = *(events_list->get(ic));
      const Trigger* trigger = event.getTrigger();
      const ASTNode* astnode = trigger->getMath();
      std::vector<std::string> dependencies_set
       = get_all_dependencies(*astnode,
                              good_params,
                              sconstant_init_assig,
                              assignment_rules_map,
                              model_reactions_map,
                              {});

      for (auto it = dependencies_set.crbegin();
        it != dependencies_set.crend(); ++it)
      {
        //genfile << "\n" << *it <<" " ;
        if (*it == "time"){
          cinitasit = sconstant_init_assig_notused.find("time_t");
          good_params.insert("time_t");
          good_params.insert("time");
          param_to_del.insert("time_t");
          param_to_del.insert("time");
        } else {
          cinitasit = sconstant_init_assig_notused.find(*it);
        }
        if (cinitasit != sconstant_init_assig_notused.cend()){
          if (cinitasit->first == "time_t"){
            wcs_const_exp.push_back("time");
            wcs_const_exp_map.insert(std::make_pair("time",cinitasit->second));
          } else {
            wcs_const_exp.push_back(cinitasit->first);
            wcs_const_exp_map.insert(std::make_pair(cinitasit->first,cinitasit->second));
          }
          good_params.insert(cinitasit->first);
          param_to_del.insert(cinitasit->first);
        }
      }
    }
  }

  // Remove parameters that have been declared
  for (const auto& x: param_to_del) {
    sconstant_init_assig.erase(x);
  }

  // Sets for keeping all consts and all vars
  std::unordered_set<std::string>::const_iterator wcs_all_const_it, wcs_all_var_it;

  // define global namespace for constants
  genfile_hdr << "namespace WCS_GLOBAL_CONST {\n";
  for (const auto& x: wcs_const) {
    genfile_hdr << "  constexpr " << Real << " " << x.first << " = " << x.second << ";\n";
    wcs_all_const.insert(x.first);
  }
  for (const auto& x: wcs_const_exp) {
    wcs_const_it = wcs_const_exp_map.find(x);
    if (wcs_const_it != wcs_const_exp_map.cend()) {
      arit = assignment_rules_map.find(SBML_formulaToString(wcs_const_it->second));
      // if there is a assignment rule for a const
      if (arit != assignment_rules_map.cend()) {
        std::set<std::string, strcomp> expr_elements;
        read_ast_node(*arit->second, expr_elements);
        unsigned int is_const = 0u;
        for (const auto& y: expr_elements) {
          wcs_const_itv = wcs_const.find(y);
          if ( wcs_const_itv != wcs_const.cend()) {
            is_const = is_const + 1;
          }
          if ( std::find(wcs_const_exp.cbegin(), wcs_const_exp.cend(), y) != wcs_const_exp.cend()) {
            is_const = is_const + 1;
          }
        }
        if (is_const == expr_elements.size()) {
          LIBSBML_CPP_NAMESPACE::ASTNode math;
          math = *arit->second;
          include_init_for_rate_rules(math, rate_rules_map);
          genfile_hdr << "  constexpr " << Real << " " << wcs_const_it->first
                  << " = " << SBML_formulaToString(&math) << ";\n";
          wcs_all_const.insert(wcs_const_it->first);
        } else {
          WCS_THROW("There is a problem with the definition of const " \
          + wcs_const_it->first + " (It takes value from a variable).");
        }

      } else {
        LIBSBML_CPP_NAMESPACE::ASTNode math;
        math = *wcs_const_it->second;
        include_init_for_rate_rules(math, rate_rules_map);
        genfile_hdr << "  constexpr " << Real << " " << wcs_const_it->first
                << " = " << SBML_formulaToString(&math) << ";\n";
        wcs_all_const.insert(wcs_const_it->first);
      }
    }
  }
  genfile_hdr << "};\n";

  // define global structure for variables
  genfile_hdr << "\n";
  genfile_hdr << "struct WCS_GLOBAL_VAR {\n";
  //insert event assignment variables
  if ( ev_assign.size() > 0ul) {
    for (const auto& x: ev_assign) {
      genfile_hdr << "  " << Real << " " << x << ";\n";
      wcs_all_var.insert(x);
    }
  }
  for (const auto& x: wcs_var) { //include events too
    wcs_all_var_it = wcs_all_var.find(x.first);
    if (wcs_all_var_it == wcs_all_var.cend()) {
      genfile_hdr << "  " << Real << " " << x.first << ";\n";
      wcs_all_var.insert(x.first);
    }
  }
  //include rate rules
  for  (const auto& x: rate_rules_map) {
    genfile_hdr << "  " << Real << " " << x.first << ";\n";
    wcs_all_var.insert(x.first);
  }
  for (const auto& x: wcs_var_exp) {
    wcs_var_it = wcs_var_exp_map.find(x);
    if (wcs_var_it != wcs_var_exp_map.cend()) {
      genfile_hdr << "  " << Real << " " << wcs_var_it->first << ";\n";
      wcs_all_var.insert(wcs_var_it->first);
    }
  }
  // print the rest of variables in assignment rules
  for (const auto& x: assignment_rules_map) {
    if ( std::find(wcs_var_exp.cbegin(), wcs_var_exp.cend(), x.first) == wcs_var_exp.cend()) {
      genfile_hdr << "  " << Real << " " << x.first <<  ";\n";
      wcs_all_var.insert(x.first);
    }
  }

  genfile_hdr << "  WCS_GLOBAL_VAR();\n";
  genfile_hdr << "  void init();\n";
  genfile_hdr << "};\n";
  genfile_hdr << "extern WCS_GLOBAL_VAR wcs_global_var;\n\n";
  genfile_impl << "WCS_GLOBAL_VAR wcs_global_var;\n\n";
  genfile_impl << "WCS_GLOBAL_VAR::WCS_GLOBAL_VAR()\n";
  genfile_impl << "{\n";
  genfile_impl << "  init();\n";
  genfile_impl << "}\n";

  genfile_impl << "void WCS_GLOBAL_VAR::init()\n";
  genfile_impl << "{\n";
  for (const auto& x: ev_assign) {
    genfile_impl << "  wcs_global_var." << x << " = WCS_GLOBAL_CONST::_init_" << x <<";\n";
  }
  for (const auto& x: wcs_var) { //include events too
    evassigit = ev_assign.find(x.first);
    if (evassigit == ev_assign.cend()) {
      genfile_impl << "  wcs_global_var." << x.first << " = " << x.second << ";\n";
    }
  }
  //include rate rules
  for  (const auto& x: rate_rules_map) {
    genfile_impl << "  wcs_global_var." << x.first << " = WCS_GLOBAL_CONST::_init_"
            << x.first << ";\n";
  }
  for (const auto& x: wcs_var_exp) {
    wcs_var_it = wcs_var_exp_map.find(x);
    if (wcs_var_it != wcs_var_exp_map.cend()) {
      LIBSBML_CPP_NAMESPACE::ASTNode math;
      math = *wcs_var_it->second;
      update_scope_ast_node(math, {}, wcs_all_const, {});
      genfile_impl << "  wcs_global_var." << wcs_var_it->first << " = "
              << SBML_formulaToString(&math) << ";\n";
    }
  }
  // print the rest of variables in assignment rules
  for (const auto& x: assignment_rules_map) {
    if ( std::find(wcs_var_exp.cbegin(), wcs_var_exp.cend(), x.first) == wcs_var_exp.cend()) {
      genfile_impl << "  wcs_global_var." << x.first <<  " = 0.0;\n";
    }
  }
  genfile_impl << "}\n";
}


void generate_cxx_code::print_functions(
  const LIBSBML_CPP_NAMESPACE::Model& model,
  std::ostream & genfile)
{
  const char* Real = generate_cxx_code::basetype_to_string<reaction_rate_t>::value;
  const ListOfFunctionDefinitions* function_definition_list
    = model.getListOfFunctionDefinitions();
  const unsigned int num_functions = function_definition_list->size();

  for (unsigned int ic = 0u; ic < num_functions; ic++) {
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
}

void generate_cxx_code::print_event_functions(
  const LIBSBML_CPP_NAMESPACE::Model& model,
  std::ostream & genfile,
  const event_assignments_t & ev_assign,
  const std::unordered_set<std::string>& wcs_all_const,
  const std::unordered_set<std::string>& wcs_all_var)
{
  const char* Real = generate_cxx_code::basetype_to_string<reaction_rate_t>::value;
  const ListOfEvents* events_list = model.getListOfEvents();
  const unsigned int num_events = events_list->size();
  for (const auto& x: ev_assign) {
    genfile << "extern \"C\" " << Real << " wcs__rate_" << x
            << "() {\n";
    genfile << "  " << Real << " " << x <<" = wcs_global_var." << x << ";\n";
    genfile << "  " << Real << " time = WCS_GLOBAL_CONST::time;\n";

    for (unsigned int ic = 0u; ic < num_events; ic++) {
      const LIBSBML_CPP_NAMESPACE::Event& event = *(events_list->get(ic));
      const ListOfEventAssignments* event_assignment_list
      = event.getListOfEventAssignments();
      const unsigned int num_event_assignments = event_assignment_list->size();
      const Trigger* trigger = event.getTrigger();
      const ASTNode* astnode = trigger->getMath();
      std::string trigger_variable = SBML_formulaToString(astnode->getChild(0));
      std::string trigger_value = SBML_formulaToString(astnode->getChild(1));
      LIBSBML_CPP_NAMESPACE::ASTNode math = *astnode->getChild(1);
      update_scope_ast_node(math, wcs_all_var, wcs_all_const, {});
      for (unsigned int ici = 0u; ici < num_event_assignments; ici++) {
        const LIBSBML_CPP_NAMESPACE::EventAssignment& eventassignment
        = *(event_assignment_list->get(ici));
        if (eventassignment.getVariable() == x){
          //genfile << "  if (variable_tr_ev == \"" << trigger_variable << "\" && value_tr_ev == "
          //        << trigger_value << ") {\n";
          genfile << "  if (time == " << SBML_formulaToString(&math) << ") {\n";
          genfile << "    " << x << " = " << SBML_formulaToString(eventassignment.getMath())
                  << ";\n";
          genfile << "    wcs_global_var." << x << " = " << x << ";\n";
          genfile << "  }\n";
        }
      }
    }

    genfile << "\n  return " << x << ";\n" << "}\n\n";
  }

}

void generate_cxx_code::print_global_state_functions(
  const LIBSBML_CPP_NAMESPACE::Model& model,
  std::ostream & genfile,
  const std::unordered_set<std::string> & good_params,
  const map_symbol_to_ast_node_t & sconstant_init_assig,
  const assignment_rules_t & assignment_rules_map,
  const model_reactions_t & model_reactions_map,
  const std::unordered_set<std::string>& wcs_all_const,
  const std::unordered_set<std::string>& wcs_all_var)
{
  const char* Real = generate_cxx_code::basetype_to_string<reaction_rate_t>::value;
  const ListOfRules* rules_list = model.getListOfRules();
  const unsigned int num_rules = rules_list->size();
  typename assignment_rules_t::const_iterator arit;
  typename model_reactions_t::const_iterator mrit;
  const std::string zero = std::string("static_cast<") + Real + ">(0)";
  for (unsigned int ic = 0u; ic < num_rules; ic++) {
    const LIBSBML_CPP_NAMESPACE::Rule& rule = *(rules_list->get(ic));
    if (rule.getType()==0) { //rate_rule
      genfile << "extern \"C\" " << Real << " wcs__rate_" << rule.getVariable()
                << "(const std::vector<" << Real << ">& __input) {\n";
      //genfile << "  int __ii=0;\n";
      std::vector<std::string> dependencies_set
        = get_all_dependencies(*rule.getMath(),
                                    good_params,
                                    sconstant_init_assig,
                                    assignment_rules_map,
                                    model_reactions_map,
                                    {});
      std::set<std::string> var_names;
      // print the function parameters with the order met in the formula
      // put first to a set all variables taking an input
      for (auto it = dependencies_set.crbegin();
           it!= dependencies_set.crend(); ++it)
      {
        arit = assignment_rules_map.find(*it);
        mrit = model_reactions_map.find(*it);
        if (arit == assignment_rules_map.cend() && mrit == model_reactions_map.cend()){
          var_names.insert(*it);
        }
      }
      std::string formula = SBML_formulaToString(rule.getMath());
      const auto input_map = build_input_map(formula, var_names);

      std::vector<std::string> var_names_ord (input_map.size());
      for (const auto& x: input_map) {
        var_names_ord[x.second] = x.first;
      }

      //print first parameters in formula taking input
      std::vector<std::string>::const_iterator it;
      int par_index = 0;
      for (it = var_names_ord.cbegin(); it < var_names_ord.cend(); it++){
        if (*it != rule.getVariable()) {
          genfile << "  " << Real << " " << *it << " = __input[" << par_index++ << "];\n";
        }
        var_names.erase(*it);
      }
      //print rest parameters (not in formula) taking input
      for (const auto& x: var_names) {
        if (x != rule.getVariable()) {
          genfile << "  " << Real << " " << x << " = __input[" << par_index++ << "];\n";
        }
      }

      LIBSBML_CPP_NAMESPACE::ASTNode math;
      // Print the rest of the parameters not taking an input
      for (auto it = dependencies_set.crbegin();
           it!= dependencies_set.crend(); ++it)
      {
        arit = assignment_rules_map.find(*it);
        mrit = model_reactions_map.find(*it);
        if (arit != assignment_rules_map.cend()){
          math = *arit->second;
          update_scope_ast_node(math, wcs_all_var, wcs_all_const, {});
          genfile << "  wcs_global_var." << arit->first << " = "
                    << SBML_formulaToString(&math) << ";\n";
        } else if (mrit != model_reactions_map.cend()){
          math = *mrit->second;
          update_scope_ast_node(math, wcs_all_var, wcs_all_const, {});
          genfile << "  " << Real << " " << mrit->first << " = "
                    << SBML_formulaToString(&math) << ";\n";
        }
      }
      math = *rule.getMath();
      update_scope_ast_node(math, wcs_all_var, wcs_all_const, {});
      genfile << "  wcs_global_var." << rule.getVariable() << " = "
              << SBML_formulaToString(&math) << ";\n";

      genfile << "  if (!isfinite(wcs_global_var." << rule.getVariable()  << ")) {\n";
      // find denominators of the dependent expressions of the rate rule formula
      for (auto it = dependencies_set.crbegin();
        it != dependencies_set.crend(); ++it)
      {
        arit = assignment_rules_map.find(*it);
        mrit = model_reactions_map.find(*it);
        if (arit != assignment_rules_map.cend()){ //assignment rules
          std::string elem_with_scope = arit->first;
          update_scope_str(elem_with_scope, wcs_all_var, wcs_all_const,{});
          genfile << "    if (!isfinite(" << elem_with_scope << ")) {\n";
          math = *arit->second;
          update_scope_ast_node(math, wcs_all_var, wcs_all_const, {});
          std::vector<std::string> denominators_noscope = return_all_denominators
          (*arit->second, rule.getVariable());
          std::vector<std::string> denominators = return_all_denominators
          (math, rule.getVariable());
          size_t sz = denominators.size();
          if (sz > 0ul) {
            for (size_t i=0ul; i<sz; i++) {
              genfile << "      if (" << denominators[i] << " == " << zero << ") {\n"
                      << "        WCS_THROW(\"Divide by zero in rate rule "
                      << rule.getVariable()  << ", in expression of " << arit->first
                      << " with the element " << denominators_noscope[i] << ".\");\n"
                      << "      }\n";
            }
          }
          genfile << "      WCS_THROW(\"Infinite or NaN result in rate rule "
                << rule.getVariable() << ", in expression " << arit->first << ".\");\n";

          genfile << "    }\n";
        } else if (mrit != model_reactions_map.cend()){ //model reactions
          std::string elem_with_scope = mrit->first;
          update_scope_str(elem_with_scope, wcs_all_var, wcs_all_const,{});
          genfile << "    if (!isfinite(" << elem_with_scope << ")) {\n";
          math = *mrit->second;
          update_scope_ast_node(math, wcs_all_var, wcs_all_const, {});
          std::vector<std::string> denominators_noscope= return_all_denominators
          (*mrit->second, rule.getVariable());
          std::vector<std::string> denominators = return_all_denominators
          (math, rule.getVariable());
          size_t sz = denominators.size();
          if (sz > 0ul) {
            for (size_t i=0ul; i<sz; i++) {
              genfile << "      if (" << denominators[i] << " == " << zero << ") {\n"
                      << "        WCS_THROW(\"Divide by zero in rate rule "
                      << rule.getVariable()  << ", in expression of " << mrit->first
                      << " with the element " << denominators_noscope[i] << ".\");\n"
                      << "      }\n";
            }
          }
          genfile << "      WCS_THROW(\"Infinite or NaN result in rate rule "
                << rule.getVariable() << ", in expression " << mrit->first << ".\");\n";

          genfile << "    }\n";
        }
      }
      // find denominators of the reaction rate formula
      math = *rule.getMath();
      update_scope_ast_node(math, wcs_all_var, wcs_all_const, {});
      std::vector<std::string> denominators_rr_noscope = return_all_denominators
      (*rule.getMath(), rule.getVariable());
      std::vector<std::string> denominators_rr = return_all_denominators
      (math, rule.getVariable());
      size_t sz = denominators_rr.size();
      if (sz > 0ul) {
        for (size_t i=0ul; i<sz; i++) {
          genfile << "    if (" << denominators_rr[i] << " == " << zero << ") {\n"
                  << "      WCS_THROW(\"Divide by zero in rate rule "
                  << rule.getVariable()
                  << " with the element " << denominators_rr_noscope[i] << ".\");\n"
                  << "    }\n";
        }
      }
      genfile << "    WCS_THROW(\"Infinite or NaN result in rate rule "
              << rule.getVariable() << ".\"); \n";
      genfile << "  }\n";




      genfile << "\n  return " << SBML_formulaToString(&math)
                << ";\n" << "}\n\n";
    }
  }
}

void generate_cxx_code::print_reaction_rates(
  const LIBSBML_CPP_NAMESPACE::Model& model,
  const unsigned int rid_start,
  const unsigned int rid_end,
  std::ostream & genfile,
  const std::string& header,
  const std::unordered_set<std::string> & good_params,
  const map_symbol_to_ast_node_t & sconstant_init_assig,
  const assignment_rules_t & assignment_rules_map,
  const model_reactions_t & model_reactions_map,
  const event_assignments_t & ev_assign,
  const std::unordered_set<std::string>& wcs_all_const,
  const std::unordered_set<std::string>& wcs_all_var,
  wcs::params_map_t& dep_params_f,
  wcs::params_map_t& dep_params_nf,
  const rate_rules_dep_t& rate_rules_dep_map)
{
  const char* Real = generate_cxx_code::basetype_to_string<reaction_rate_t>::value;
  const std::string zero = std::string("static_cast<") + Real + ">(0)";
  const ListOfReactions* reaction_list = model.getListOfReactions();
  const unsigned int num_reactions = reaction_list->size();
  typename assignment_rules_t::const_iterator arit;
  wcs::rate_rules_dep_t::const_iterator rrdit;
  using  event_assignments_t = std::unordered_set<std::string>;
  typename event_assignments_t::const_iterator evassigit;

  if (rid_start > num_reactions ||
      rid_end > num_reactions ||
      rid_start > rid_end) {
    WCS_THROW("Invalid reaction range: " + std::to_string(rid_start) + \
              " - "+ std::to_string(rid_end));
  }

  genfile << "#include \"" + header + '"' + "\n\n";
  for (unsigned int ic = rid_start; ic < rid_end; ic++) {
    const LIBSBML_CPP_NAMESPACE::Reaction& reaction = *(reaction_list->get(ic));
    // const LIBSBML_CPP_NAMESPACE::ListOfLocalParameters* local_parameter_list
    //   = reaction.getKineticLaw()->getListOfLocalParameters();
    // unsigned int num_localparameters = local_parameter_list->size();
    // a map to keep the compartment and their units
    using  compartment_t = std::unordered_map <std::string,double>;
    typename compartment_t::const_iterator rci;
    compartment_t reaction_comp;
    unsigned int num_localparameters = 0u;
    if (model.getLevel() > 2) {
      const LIBSBML_CPP_NAMESPACE::ListOfLocalParameters* local_parameter_list
        = reaction.getKineticLaw()->getListOfLocalParameters();
      num_localparameters = local_parameter_list->size();
    } else {
      const LIBSBML_CPP_NAMESPACE::ListOfParameters* local_parameter_list
        = reaction.getKineticLaw()->getListOfParameters();
      num_localparameters = local_parameter_list->size(); 
    }  

    
    genfile << "extern \"C\" " << Real << " wcs__rate_" << reaction.getIdAttribute()
              << "(const std::vector<" << Real <<">& __input) {\n";
              
    double sum_stoich = 0.0;
    //print reaction's local parameters
    using reaction_local_parameters_t = std::unordered_set<std::string>;
    reaction_local_parameters_t localpset;
    // for (unsigned int pi = 0u; pi < num_localparameters; pi++) {
    //   const LIBSBML_CPP_NAMESPACE::LocalParameter* localparameter = local_parameter_list->get(pi);
    //   genfile << "  constexpr " << Real << " " << localparameter->getIdAttribute() << " = "
    //           << localparameter->getValue() << ";\n";
    //   localpset.insert(localparameter->getIdAttribute());
    if (model.getLevel() > 2) {
      const LIBSBML_CPP_NAMESPACE::ListOfLocalParameters* local_parameter_list
        = reaction.getKineticLaw()->getListOfLocalParameters();
      for (unsigned int pi = 0u; pi < num_localparameters; pi++) {
        const LIBSBML_CPP_NAMESPACE::LocalParameter* localparameter = local_parameter_list->get(pi);
        genfile << "  constexpr " << Real << " " << localparameter->getIdAttribute() << " = "
               << localparameter->getValue() << ";\n";
        localpset.insert(localparameter->getIdAttribute());
      }
    } else {
      const LIBSBML_CPP_NAMESPACE::ListOfParameters* local_parameter_list
        = reaction.getKineticLaw()->getListOfParameters();
      for (unsigned int pi = 0u; pi < num_localparameters; pi++) {
        const LIBSBML_CPP_NAMESPACE::Parameter* localparameter = local_parameter_list->get(pi);
        genfile << "  constexpr " << Real << " " << localparameter->getIdAttribute() << " = "
               << localparameter->getValue() << ";\n";
        localpset.insert(localparameter->getIdAttribute());
      }
    }
    //genfile << "  int __ii=0;\n";
    //genfile << "printf(\"Number of input arguments of %s received : %lu \", \""
    //        << reaction.getIdAttribute() << "\", __input.size());\n";
    //genfile << "printf(\"nfunction: %s \\n\", \"" << reaction.getIdAttribute() << "\");\n";
    std::vector<std::string> dependencies_set
    = get_all_dependencies(*reaction.getKineticLaw()->getMath(),
                          good_params,
                          sconstant_init_assig,
                          assignment_rules_map,
                          model_reactions_map,
                          localpset);

    std::set<std::string> var_names, all_var_names;
    std::set<std::string>::const_iterator var_names_it, all_var_names_it;

    // print the function parameters with the order met in the formula
    //put first to a set all variables taking an input
    for (auto it = dependencies_set.crbegin();
         it != dependencies_set.crend(); ++it)
    {
      arit = assignment_rules_map.find(*it);
      if (arit == assignment_rules_map.cend()){
        var_names.insert(*it);
      }
    }
    std::string formula = SBML_formulaToString(reaction.getKineticLaw()->getMath());
    const auto input_map = build_input_map(formula, var_names);

    std::vector<std::string> var_names_ord (input_map.size());
    for (const auto& x: input_map) {
      var_names_ord[x.second] = x.first;
    }
    //print first parameters in formula taking input
    std::vector<std::string>::const_iterator it;
    std::vector<std::string> par_names, par_names_nf;
    all_var_names = var_names;
    int par_index = 0;
    const ListOfSpecies* species_list = model.getListOfSpecies();
    double avog_num = 6.02214e+23;
    bool concentration_used = false;

    for (it = var_names_ord.cbegin(); it < var_names_ord.cend(); it++){
      rrdit = rate_rules_dep_map.find(*it);
      evassigit = ev_assign.find(*it);
      if (rrdit != rate_rules_dep_map.cend()) { //rate rules
        const std::set<std::string>& params_fn = rrdit->second;
        std::set<std::string> params_fn_dec;
        std::string function_input;
        std::set<std::string>::const_iterator itf;
        // add all the dependent species of rate rules
        for (itf = params_fn.cbegin(); itf != params_fn.cend(); itf++){
          if (itf != params_fn.cbegin()) {
           function_input = function_input + ", ";
          }
          var_names_it = var_names.find(*itf);
          all_var_names_it = all_var_names.find(*itf);
          if (all_var_names_it == all_var_names.cend()) {
            genfile << "  " << Real << " " << *itf << " = __input[" << par_index++ << "];\n";
            if (species_list->get(*itf) != NULL) {
              if (species_list->get(*itf)->isSetInitialConcentration()) {
                double comp_unit;
                std::string compart_name; 
                double convert_factor = count2Concentration_converter(model, *species_list->get(*itf), compart_name, comp_unit);
                // Add values for the compartment in map
                rci = reaction_comp.find(compart_name) ;
                if (rci == reaction_comp.end()) {
                  reaction_comp.insert(std::make_pair(compart_name, comp_unit));
                }
                genfile << "  " << *itf << " = " <<  *itf << " / " << convert_factor <<  ";\n";
                concentration_used = true;
              }
            }
            var_names.erase(*itf);
          } else {
            if (var_names_it != var_names.cend()) {
              genfile << "  " << Real << " " << *itf << " = __input[" << par_index++ << "];\n";
              if (species_list->get(*itf) != NULL) {
                if (species_list->get(*itf)->isSetInitialConcentration()) {
                  double comp_unit;
                  std::string compart_name; 
                  double convert_factor = count2Concentration_converter(model, *species_list->get(*itf), compart_name, comp_unit);
                  // Add values for the compartment in map
                  rci = reaction_comp.find(compart_name) ;
                  if (rci == reaction_comp.end()) {
                    reaction_comp.insert(std::make_pair(compart_name, comp_unit));
                  }
                  genfile << "  " << *itf << " = " <<  *itf << " / " << convert_factor <<  ";\n";
                  concentration_used = true;
                }
              }
              var_names.erase(*itf);
            }
          }
          function_input = function_input + *itf ;
        }
        genfile << "  wcs_global_var."  << *it << " = wcs__rate_" << *it << "(std::vector<"
                << Real << "> {" << function_input << "});\n";
        par_names.push_back(*it);
      } else if (evassigit != ev_assign.cend()) { //events variables
        std::string function_input;
        function_input ="";
        genfile << "  wcs_global_var." << *it << " = wcs__rate_" << *it << "(" << function_input << ");\n";
      } else {
        genfile << "  " << Real << " " << *it << " = __input[" << par_index++ << "];\n";
        if (species_list->get(*it) != NULL) {
          if (species_list->get(*it)->isSetInitialConcentration()) {
            double comp_unit;
            std::string compart_name; 
            double convert_factor = count2Concentration_converter(model, *species_list->get(*it), compart_name, comp_unit);
            // Add values for the compartment in map
            rci = reaction_comp.find(compart_name) ;
            if (rci == reaction_comp.end()) {
              reaction_comp.insert(std::make_pair(compart_name, comp_unit));
            }
            genfile << "  " << *it << " = " <<  *it << " / " << convert_factor <<  ";\n";
            concentration_used = true;
          }
        }
        if (reaction.getReactant(*it) != NULL) {
          sum_stoich = sum_stoich + reaction.getReactant(*it)->getStoichiometry();
        }
        par_names.push_back(*it);
      }
      var_names.erase(*it);
    }
    //print rest parameters (not in formula) taking input
    for (const auto& x: var_names) {
      rrdit = rate_rules_dep_map.find(x);
      evassigit = ev_assign.find(x);
      if (rrdit != rate_rules_dep_map.cend()) { //rate rules
        const std::set<std::string>& params_fn = rrdit->second;
        std::string function_input ;
        std::set<std::string>::const_iterator itf;
        // add all the dependent species of rate rules
        for (itf = params_fn.cbegin(); itf != params_fn.cend(); itf++){
          if (itf != params_fn.cbegin()) {
            function_input = function_input + ", ";
          }
          var_names_it = var_names.find(*itf);
          all_var_names_it = all_var_names.find(*itf);
          if (all_var_names_it == all_var_names.cend()) {
              genfile << "  " << Real << " " << *itf << " = __input[" << par_index++ << "];\n";
              var_names.erase(*itf);
          } else {
            if (var_names_it != var_names.cend()) {
              genfile << "  " << Real << " " << *itf << " = __input[" << par_index++ << "];\n";
              var_names.erase(*itf);
            }
          }
          function_input = function_input + *itf ;
        }
        genfile << "  wcs_global_var."  << x << " = wcs__rate_" << x << "(std::vector<" << Real
                << "> {" << function_input << "});\n";
        par_names_nf.push_back(x);
      } else if (evassigit != ev_assign.cend()) { // event variables
        std::string function_input;
        function_input ="";
        genfile << "  wcs_global_var." << x << " = wcs__rate_" << x << "(" << function_input << ");\n";
      } else {
        genfile << "  " << Real << " " << x << " = __input[" << par_index++ << "];\n";
        par_names_nf.push_back(x);
      }
    }
    // put input parameters into maps
    dep_params_f.insert(std::make_pair(reaction.getIdAttribute(),
                                       par_names));
    dep_params_nf.insert(std::make_pair(reaction.getIdAttribute(),
                                        par_names_nf));

    //genfile << "printf(\" and expected  %u \\n\", " << par_index << ");\n";
    /*genfile << "  printf(\"Expected in generated code: \");\n";
    for (auto& x: par_names) {
      genfile << "  printf(\"%s, \", \"" << x << "\");\n";
    }
    for (auto& x: par_names_nf) {
      genfile << "  printf(\" $ %s, \", \"" << x << "\");\n";
    }*/

    std::set<LIBSBML_CPP_NAMESPACE::ASTNode> all_denominators;
    LIBSBML_CPP_NAMESPACE::ASTNode math;

    //print the rest of the parameters not taking an input
    for (auto it = dependencies_set.crbegin();
         it != dependencies_set.crend(); ++it)
    {
      arit = assignment_rules_map.find(*it);
      if (arit != assignment_rules_map.cend()){
        math = *arit->second;
        update_scope_ast_node(math, wcs_all_var, wcs_all_const, localpset);
        genfile << "  wcs_global_var." << arit->first << " = "
                  << SBML_formulaToString( &math) << ";\n";
      }
    }
    math = *reaction.getKineticLaw()->getMath();
    update_scope_ast_node(math, wcs_all_var, wcs_all_const, localpset);
    std::string new_math = sbml_functions::SBML_formulaToString_wcs(&math);
    genfile << "  " << Real << " " << reaction.getIdAttribute() << " = "
              << new_math //SBML_formulaToString( &math)
              << ";\n";
    genfile << "  if (!isfinite(" << reaction.getIdAttribute() << ")) {\n";
    // find denominators of the dependent expressions of the reaction rate formula
    for (auto it = dependencies_set.crbegin();
         it != dependencies_set.crend(); ++it)
    {
      arit = assignment_rules_map.find(*it);
      if (arit != assignment_rules_map.cend()){
        std::string elem_with_scope = arit->first;
        update_scope_str(elem_with_scope, wcs_all_var, wcs_all_const,localpset);
        genfile << "    if (!isfinite(" << elem_with_scope << ")) {\n";
        math = *arit->second;
        update_scope_ast_node(math, wcs_all_var, wcs_all_const, localpset);
        std::vector<std::string> denominators_noscope =
        return_all_denominators (*arit->second, reaction.getIdAttribute());
        std::vector<std::string> denominators = return_all_denominators
        (math, reaction.getIdAttribute());
        size_t sz = denominators.size();
        if (sz > 0ul) {
          for (size_t i=0ul; i<sz; i++) {
            genfile << "      if (" << denominators[i] << " == " << zero << ") {\n"
                    << "        WCS_THROW(\"Divide by zero in reaction "
                    << reaction.getIdAttribute() << ", in expression of " << arit->first
                    << " with the element " << denominators_noscope[i] << ".\"); \n"
                    << "      }\n";
          }
        }
        genfile << "      WCS_THROW(\"Infinite or NaN result in reaction "
                << reaction.getIdAttribute() << ", in expression " << arit->first << ".\"); \n";

        genfile << "    }\n";
      }
    }
    // find denominators of the reaction rate formula
    math = *reaction.getKineticLaw()->getMath();
    update_scope_ast_node(math, wcs_all_var, wcs_all_const, localpset);
    std::vector<std::string> denominators_rr_noscope = return_all_denominators
    (*reaction.getKineticLaw()->getMath(), reaction.getIdAttribute());
    std::vector<std::string> denominators_rr = return_all_denominators
    (math, reaction.getIdAttribute());
    size_t sz = denominators_rr.size();
    if (sz > 0ul) {
      for (size_t i=0ul; i<sz; i++) {
        genfile << "    if (" << denominators_rr[i] << " == " << zero << ") {\n"
                << "      WCS_THROW(\"Divide by zero in reaction "
                << reaction.getIdAttribute()
                << " with the element " << denominators_rr_noscope[i] << ".\"); \n"
                << "    }\n";
      }
    }
    genfile << "    WCS_THROW(\"Infinite or NaN result in reaction "
            << reaction.getIdAttribute() << ".\"); \n";
    genfile << "  }\n";
    //genfile << "  printf(\" Result of %s : %f \", \"" << reaction.getIdAttribute() << "\","
    //        << reaction.getIdAttribute() << ");\n";
    if (!concentration_used) {
      genfile << "  return " << reaction.getIdAttribute() << ";\n";
    } else {
      // Calculate the compartment unit for all the compartments
      double comp_unit = 1.0;
      for ( auto rci = reaction_comp.begin(); rci != reaction_comp.end(); ++rci ){
        comp_unit = comp_unit * rci->second; 
      }
      genfile << "  return " << reaction.getIdAttribute() << " * " << avog_num << " * " << comp_unit << ";\n"; 
    }
    genfile << "}\n\n";
  }
}


/**
 *  If `regen` is set to false (which is the default), then the library file
 *  at the given path is reused. If no file exists at the path specified, or
 *  the extention of the file is not `.so`, which is for a dynamic library,
 *  then a new file is generated with the same name except that the extension
 *  is replaced with `.so`. If the given path is an empty string, then
 *  `[wcs_default_gen_name].so` file is generated under current directory.
 *  Caller can set `save_log` for dignosis if compilation failure is expected.
 *  By default, it is off. Caller can also set or unset `cleanup` argument
 *  to remove the temporary source file generated or leave it.
 *  By default, it is on. To avoid generating a single translation unit that
 *  is extreme large, we limit the number of reactions that are written
 *  into a single file by specifying chunk_size. By default, it is set to
 *  1000.
 */
generate_cxx_code::generate_cxx_code(const std::string& libpath,
                                     bool regen, bool save_log, bool cleanup,
                                     const std::string& tmp_dir,
                                     unsigned int chunk_size,
                                     unsigned int num_compiling_threads)
: m_lib_filename(libpath), m_regen(regen), m_save_log(save_log),
  m_cleanup(cleanup), m_tmp_dir(tmp_dir), m_chunk(chunk_size),
  m_num_compiling_threads(num_compiling_threads)
{
  m_regen = m_regen || !check_if_file_exists(m_lib_filename);

  if (m_lib_filename.empty()) {
    m_regen = true;
    m_lib_filename = std::string(wcs_default_gen_name) + ".so";
  }

  std::string dir;
  std::string stem;
  std::string ext;

  extract_file_component(m_lib_filename, dir, stem, ext);

  if (ext != ".so") {
    m_regen = true;
    m_lib_filename = dir + stem + ".so";
  }
  if (m_lib_filename[0] != '/') {
     char canonical[PATH_MAX] = {'\0'};
     realpath("./", canonical);
     const size_t sz = strlen(canonical);
     if ((sz > 0ul) && (canonical[sz-1] != '/')) {
       m_lib_filename = canonical + std::string{"/"} + m_lib_filename;
     } else { // it is ok even if canonical is empty
       m_lib_filename = canonical + m_lib_filename;
     }
  }

 #if defined(_OPENMP)
  #pragma omp master
 #endif // defined(_OPENMP)
  {
    if (((m_tmp_dir.size() == 1u) && (m_tmp_dir[0] == '/')) ||
        m_tmp_dir.empty())
    {
      std::cerr << "Creating files into the current directory" << std::endl;
      m_tmp_dir = ".";
    } else if (m_tmp_dir != "/tmp") {
      if (m_regen) {
        int ret = mkdir_as_needed(m_tmp_dir);
        if (ret != 0) {
          WCS_THROW("Terminating after failed to create a directory.");
        }
      }
    }
  }
 #if defined(_OPENMP)
  #pragma omp barrier
 #endif // defined(_OPENMP)

  if (m_tmp_dir[0] != '/') {
     char canonical[PATH_MAX] = {'\0'};
     realpath("./", canonical);
     const size_t sz = strlen(canonical);
     if ((sz > 0ul) && (canonical[sz-1] != '/')) {
       m_tmp_dir = canonical + std::string{"/"} + m_tmp_dir;
     } else { // it is ok even if canonical is empty
       m_tmp_dir = canonical + m_tmp_dir;
     }
  }

 #if defined(_OPENMP)
  #pragma omp master
 #endif // defined(_OPENMP)
  {
    std::cerr << std::string(m_regen? "Generating" : "Reusing")
              + " the machine code for reaction rate formula ("
              + m_lib_filename + ")" << std::endl;
  }
}

//https://stackoverflow.com/questions/8243743/is-there-a-null-stdostream-implementation-in-c-or-libraries
class nullstream : public std::ostream {
 public:
  nullstream() : std::ostream(nullptr) {}
};

template <typename T>
const nullstream &operator<<(nullstream &&os, const T&) {
  return os;
}

void generate_cxx_code::create_ostream(src_file_t& ofile, size_t suffix_size)
{
  std::string& src_name = ofile.first;
  std::unique_ptr<std::ostream>& os_ptr = ofile.second;


  char tmp_filename[PATH_MAX+1] = {'\0'};
  strncpy(tmp_filename, src_name.c_str(), PATH_MAX);

  int fd = mkostemps(tmp_filename, static_cast<int>(suffix_size),
                     O_CREAT | O_RDWR | O_EXCL);
  std::stringstream ss;
  ss << tmp_filename;
  src_name = ss.str();

  if (fd < 1) {
    WCS_THROW("\n Creation of temp file '" + src_name +
              "' failed with error " + strerror(errno));
  }

  close(fd);

  os_ptr = std::make_unique<std::ofstream>(src_name);
  if (!os_ptr || !(*os_ptr)) {
    WCS_THROW("\n Failed to open a source file " + src_name);
  }
}

/**
 *  Open an output stream for the code to be generated. In case that, it is
 *  set to generate code, open a temporary file with a unique name under /tmp,
 *  which has an extension `.cc`. If it is not set to generate a code, but
 *  to reuse the existing library file, a null stream is open.
 */
void generate_cxx_code::open_ostream(const unsigned int num_reactions)
{
 #if (__GLIBC__ < 2) || (__GLIBC_MINOR__ < 11)
   #error Requires Glibc version >= 2.19 for mkostemps()
   //std::cerr << "Glibc version: " << __GLIBC__ << "." << __GLIBC_MINOR__ << std::endl;
 #endif

  if (m_chunk == 0u) {
    m_chunk = 3000u;
  }
  size_t num_reaction_files = (num_reactions + m_chunk - 1) / m_chunk;
  m_ostreams.resize(num_reaction_files + 3);

  for (auto& os: m_ostreams) {
    os.first = "";
    os.second = std::make_unique<nullstream>();
  }

  if (m_regen) {
   /* In case of running OpenMP with partitioned network, each thread has
    * its own generate_cxx_code object, we need to prevent redundant file
    * I/O. However, each thread still requires dependency analysis.
    */
   #if defined(_OPENMP)
    m_regen = false; // default for non-master
    #pragma omp barrier
    #pragma omp master
   #endif // defined(_OPENMP)
    { // only the master executes this block in case of parallel execution
      m_regen = true; // only master opens/closes a stream
      std::string dir;
      std::string stem;
      std::string ext;
      extract_file_component(m_lib_filename, dir, stem, ext);

      m_ostreams.clear();
      m_ostreams.resize(num_reaction_files + 3);

      const std::string hdr_suffix = ".hpp";
      const std::string src_suffix = ".cpp";
      m_ostreams[0].first = m_tmp_dir + "/Makefile_XXXXXX";
      m_ostreams[1].first = m_tmp_dir + "/" + stem + "_XXXXXX" + hdr_suffix;
      m_ostreams[2].first = m_tmp_dir + "/" + stem + "_XXXXXX" + src_suffix;
      create_ostream(m_ostreams[0], 0);
      create_ostream(m_ostreams[1], hdr_suffix.length());
      create_ostream(m_ostreams[2], src_suffix.length());

      for (unsigned i = 0u, j = 0u; i < num_reactions; i += m_chunk, j++) {
        m_ostreams[j+3].first = m_tmp_dir + "/" + stem + '_' + std::to_string(j)
                              + "_XXXXXX" + src_suffix;
        create_ostream(m_ostreams[j+3], src_suffix.size());
      }
    }
   #if defined(_OPENMP)
    #pragma omp barrier
   #endif // defined(_OPENMP)
  }
}

void generate_cxx_code::close_ostream(std::unique_ptr<std::ostream>& os_ptr)
{
  if (m_regen) {
   #if defined(_OPENMP)
    // avoid closing the file multiple times
    #pragma omp master
   #endif // defined(_OPENMP)
    {
      if (dynamic_cast<std::ofstream*>(os_ptr.get()) != nullptr) {
        fsync_ofstream(dynamic_cast<std::ofstream&>(*os_ptr));
        dynamic_cast<std::ofstream*>(os_ptr.get())->flush();
        dynamic_cast<std::ofstream*>(os_ptr.get())->close();
        delete os_ptr.release();
        os_ptr = nullptr;
      }
    }
  }
}

static unsigned int get_model_species(
  const LIBSBML_CPP_NAMESPACE::Model& model,
  assignment_rules_t& assignment_rules_map,
  std::unordered_set <std::string>& model_species)
{
  const ListOfSpecies* species_list = model.getListOfSpecies();
  const unsigned int num_species = species_list->size();

  for (unsigned int ic = 0u; ic < num_species; ic++) {
    if (assignment_rules_map.find(species_list->get(ic)->getIdAttribute())
        == assignment_rules_map.cend())
    {
      model_species.insert(species_list->get(ic)->getIdAttribute());
    }
  }
  return num_species;
}

static unsigned int get_sinitial_assignments(
  const LIBSBML_CPP_NAMESPACE::Model& model,
  initial_assignments_t& sinitial_assignments)
{
  const ListOfInitialAssignments* assignments_list
    = model.getListOfInitialAssignments();
  const unsigned int num_assignments = assignments_list->size();

  for (unsigned int ic = 0u; ic < num_assignments; ic++) {
    const LIBSBML_CPP_NAMESPACE::InitialAssignment& initialassignment
      = *(assignments_list->get(ic));
    sinitial_assignments.insert(std::make_pair(initialassignment.getSymbol(),
                                               initialassignment.getMath()));
  }
  return num_assignments;
}

static unsigned int get_model_reactions_map(
  const LIBSBML_CPP_NAMESPACE::Model& model,
  model_reactions_t& model_reactions_map)
{
  const ListOfReactions* reaction_list = model.getListOfReactions();
  const unsigned int num_reactions = reaction_list->size();

  for (unsigned int ic = 0; ic < num_reactions; ic++) {
    const LIBSBML_CPP_NAMESPACE::Reaction& reaction = *(reaction_list->get(ic));

    if (!reaction.isSetKineticLaw()) {
      WCS_THROW("The formula of the reaction " + reaction.getIdAttribute() \
                + " should be set.");
    }

    model_reactions_map.insert(std::make_pair(reaction.getIdAttribute(),
                                              reaction.getKineticLaw()->getMath()));
  }
  return num_reactions;
}

static unsigned int get_ev_assign(
  const LIBSBML_CPP_NAMESPACE::Model& model,
  event_assignments_t& ev_assign)
{
  const ListOfEvents* events_list = model.getListOfEvents();
  const unsigned int num_events = events_list->size();

  for (unsigned int ic = 0u; ic < num_events; ic++) {
    const LIBSBML_CPP_NAMESPACE::Event& event = *(events_list->get(ic));
    const ListOfEventAssignments* event_assignment_list
      = event.getListOfEventAssignments();
    const unsigned int num_event_assignments = event_assignment_list->size();

    for (unsigned int ici = 0u; ici < num_event_assignments; ici++) {
      const LIBSBML_CPP_NAMESPACE::EventAssignment& eventassignment
        = *(event_assignment_list->get(ici));
      if (ev_assign.find(eventassignment.getVariable()) == ev_assign.cend()) {
        //genfile << eventassignment.getVariable() <<"\n";
        ev_assign.insert(eventassignment.getVariable());
      }
    }
  }
  return num_events;
}

static unsigned int get_rules_map(
  const LIBSBML_CPP_NAMESPACE::Model& model,
  rate_rules_t rate_rules_map,
  assignment_rules_t& assignment_rules_map)
{
  const ListOfRules* rules_list = model.getListOfRules();
  const unsigned int num_rules = rules_list->size();

  for (unsigned int ic = 0; ic < num_rules; ic++) {
    const LIBSBML_CPP_NAMESPACE::Rule& rule = *(rules_list->get(ic));
    if (rule.getType() == 0) { //rate_rule
      rate_rules_map.insert(std::make_pair(rule.getVariable(), rule.getMath()));
    }
    if (rule.getType() == 1) { //assignment_rule
      assignment_rules_map.insert(std::make_pair(rule.getVariable(), rule.getMath()));
    }
  }
  return num_rules;
}

void generate_cxx_code::get_rate_rules_dep_map(
  const rate_rules_t& rate_rules_map,
  const assignment_rules_t& assignment_rules_map,
  const std::unordered_set<std::string>& good_params,
  const constant_init_ass_t& sconstant_init_assig,
  const model_reactions_t& model_reactions_map,
  rate_rules_dep_t& rate_rules_dep_map)
{
  for  (const auto& x: rate_rules_map) {
    std::set<std::string> var_f;

    std::vector<std::string> dependencies_set
      = generate_cxx_code::get_all_dependencies(
                             *x.second,
                             good_params,
                             sconstant_init_assig,
                             assignment_rules_map,
                             model_reactions_map,
                             {});

    for (auto it = dependencies_set.crbegin();
         it!= dependencies_set.crend(); ++it)
    {
      typename assignment_rules_t::const_iterator arit
        = assignment_rules_map.find(*it);
      typename model_reactions_t::const_iterator mrit
        = model_reactions_map.find(*it);

      if (arit == assignment_rules_map.cend() && mrit == model_reactions_map.cend())
      {
        if  (*it != x.first) {
          var_f.insert(*it);
        }
      }
    }

    rate_rules_dep_map.insert(std::make_pair(x.first, var_f));
  }
}

void generate_cxx_code::write_header(
  const LIBSBML_CPP_NAMESPACE::Model& model,
  std::ostream& genfile)
{
  const ListOfFunctionDefinitions* function_definition_list
    = model.getListOfFunctionDefinitions();
  const unsigned int num_functions = function_definition_list->size();

  const char * Real = generate_cxx_code::basetype_to_string<reaction_rate_t>::value;

  genfile << "/** Autogenerated source code, do not edit! */\n"
            << "#ifndef WCS_REACTION_RATE_EVALUTATION_FUNCTIONS_GENERATED\n"
            << "#define WCS_REACTION_RATE_EVALUTATION_FUNCTIONS_GENERATED\n"
            << "\n//C++ includes\n"
            << "#include <vector>\n"
            << "#include <cmath>\n"
            << "#include <cstdio>\n"
            << "#include <string>\n"
            << "#include <math.h>\n"
            << "#include <iostream>\n"
            << "#include <algorithm>\n"
            << "#include \"utils/exception.hpp\"\n\n"
            << "//Include the text of the SBML in here.\n"
            << "//Use \"xxd -i\" to get the text as a C symbol.\n"
            << "extern \"C\" \n{\n"
            << "  extern unsigned char __original_sbml[];\n}\n\n"
            << "//Constants and other functions defined by the SBML standard.\n\n"
            << "//Get the correct floating point type from the code at runtime.\n"
            << "typedef " << Real << " reaction_rate_t;\n"
            //<< "typedef reaction_rate_t " << Real << ";\n\n"
            << "//Prototype all the functions\n";

  for (unsigned int ic = 0u; ic < num_functions; ic++) {
    const LIBSBML_CPP_NAMESPACE::FunctionDefinition& function
      = *(function_definition_list->get(ic));
    genfile << "static inline " << Real << " " << function.getIdAttribute() <<"(";
    unsigned int num_arg = function.getNumArguments();

    for (unsigned int ia = 0u; ia < num_arg; ia++) {
      genfile << Real << " " << function.getArgument(ia)->getName();
      if (ia!=num_arg-1) {
        genfile << ", ";
      }
    }
    genfile << ");\n";
  }
  genfile << "\n";
  genfile << "//Define all the constants and initial states\n";
}

void generate_cxx_code::write_common_impl(
  const LIBSBML_CPP_NAMESPACE::Model& model,
  const std::string& header_name,
  std::ostream& genfile)
{
  const ListOfFunctionDefinitions* function_definition_list
    = model.getListOfFunctionDefinitions();
  const unsigned int num_functions = function_definition_list->size();

  genfile << "#include \"" + header_name + '"' + "\n\n";
  genfile << "extern \"C\" \n{\n"
          << "  unsigned char __original_sbml[] = \"\";\n}\n\n";

  if ( num_functions > 0ul) {
    genfile << "\n//Define the functions\n";
    generate_cxx_code::print_functions(model, genfile);
  }
}


void generate_cxx_code::generate_code(
  const LIBSBML_CPP_NAMESPACE::Model& model,
  params_map_t& dep_params_f,
  params_map_t& dep_params_nf,
  rate_rules_dep_t& rate_rules_dep_map)
{
 #if defined(_OPENMP)
  #pragma omp master
 #endif // defined(_OPENMP)
  {
    std::string msg = std::string("Analyzing dependencies ")
                    + (m_regen? " and generating code for JIT ..." : "...");
    std::cerr << msg << std::endl;
  }

  // A map for initial_assignments
  initial_assignments_t sinitial_assignments;

  // Put initial assignments in map
  get_sinitial_assignments(model, sinitial_assignments);

  // A map for model reactions
  model_reactions_t model_reactions_map;

  // Put reactions in a map
  const unsigned int num_reactions
    = get_model_reactions_map(model, model_reactions_map);

  // A set for model event assignments
  event_assignments_t ev_assign;

  //Put event assignments in a set
  get_ev_assign(model, ev_assign);

  // A map for model rate rules
  rate_rules_t rate_rules_map;

  // A map for model assignment rules
  assignment_rules_t assignment_rules_map;

  // Put rate rules and assignements rules in maps
  get_rules_map(model, rate_rules_map, assignment_rules_map);

  std::unordered_set <std::string> model_species;

  // Put species in a set
  get_model_species(model, assignment_rules_map, model_species);
  std::unordered_set <std::string> used_params;

  // Find used parameters in the rates and the differential rates
  // The rest aruments are read-only
  generate_cxx_code::find_used_params(
    model, used_params, sinitial_assignments, assignment_rules_map,
    model_reactions_map, rate_rules_map, model_species, ev_assign);


  open_ostream(num_reactions);
  std::ostream& os_header = *(m_ostreams[1].second);
  std::ostream& os_common_impl = *(m_ostreams[2].second);


  write_header(model, os_header);
  write_common_impl(model, m_ostreams[1].first, os_common_impl);

  //  A map for constants in initial assignments
  constant_init_ass_t sconstant_init_assig;
  std::unordered_set<std::string> good_params;
  std::unordered_set<std::string> wcs_all_const, wcs_all_var;

  // Populate sconstant_init_assig, good_params, wcs_all_const, and wcs_all_var
  generate_cxx_code::print_constants_and_initial_states(
    model, os_header, os_common_impl,
    sconstant_init_assig, sinitial_assignments,
    assignment_rules_map, used_params, good_params, model_reactions_map,
    rate_rules_map, wcs_all_const, wcs_all_var, ev_assign);


  // Put dependencies of rate rules in a map (for transient parameters)
  get_rate_rules_dep_map(rate_rules_map, assignment_rules_map, good_params,
                         sconstant_init_assig, model_reactions_map,
                         rate_rules_dep_map);

  // define functions for events
  if (ev_assign.size() > 0ul) {
    os_common_impl << "\n//Define functions for events\n";
    generate_cxx_code::print_event_functions(model, os_common_impl, ev_assign,
                                             wcs_all_const, wcs_all_var);
  }

  os_header << "\n//Define the functions for updating global state variables\n";
  os_common_impl << "\n//Define the functions for updating global state variables\n";
  generate_cxx_code::print_global_state_functions(
    model, os_common_impl, good_params, sconstant_init_assig,
    assignment_rules_map, model_reactions_map, wcs_all_const, wcs_all_var);

  os_header << "\n#endif // WCS_REACTION_RATE_EVALUTATION_FUNCTIONS_GENERATED\n";
  close_ostream(m_ostreams[1].second);
  close_ostream(m_ostreams[2].second);

  for (unsigned i = 0u, j = 0u; i < num_reactions; i += m_chunk, j++) {
    std::ostream& genfile = *(m_ostreams[j+3].second);
    genfile << "//Define the rates\n";
    const unsigned int rid_end = std::min(i + m_chunk, num_reactions);

    generate_cxx_code::print_reaction_rates(
      model, i, rid_end, genfile, m_ostreams[1].first,
      good_params, sconstant_init_assig,
      assignment_rules_map, model_reactions_map, ev_assign,
      wcs_all_const, wcs_all_var, dep_params_f,
      dep_params_nf, rate_rules_dep_map);
    close_ostream(m_ostreams[j+3].second);
  }
}

static int build(const std::string& cmd,
                 const std::string& target_name,
                 const std::string& tmp_file,
                 const std::string& compilation_log)
{
  int ret = system(cmd.c_str());
  if (!WIFEXITED(ret)) {
    std::string msg = "The command `" + cmd + "` has failed.";
    std::cerr << msg << std::endl;
    ret = EXIT_FAILURE;
  } else if (WEXITSTATUS(ret) != 0) {
    std::string msg = "Failed to build " + target_name + ".";
                    + " The command used is\n" + cmd;
    if (!compilation_log.empty()) {
      msg += " See wcs_jit_log.txt for further details.";
    }
    std::cerr << msg << std::endl;
    ret = EXIT_FAILURE;
  }

  if (!tmp_file.empty()) {
    std::string cmd_rm = "rm -f " + tmp_file;
    ret = system(cmd_rm.c_str());
    if (!WIFEXITED(ret)) {
      std::string msg = "The command `" + cmd_rm + "` has failed.";
      std::cerr << msg << std::endl;
      ret = EXIT_FAILURE;
    }
  }
  return ret;
}

std::string generate_cxx_code::gen_makefile()
{
  std::string obj_files;

 #if defined(_OPENMP)
  #pragma omp master
 #endif // defined(_OPENMP)
  {
    const std::string& hdr_filename
      = get_subpath(m_tmp_dir, m_ostreams[1].first);
    std::ostream& os_makefile = *(m_ostreams[0].second);
    std::string shared_lib;

    {
      std::string dir, stem, ext;
      extract_file_component(m_lib_filename, dir, stem, ext);
      shared_lib = stem + ext;

      os_makefile << "CXX = " + std::string(CMAKE_CXX_COMPILER) + "\n";
      os_makefile << "CXXFLAGS = " + std::string(CMAKE_CXX_FLAGS) + "\n";
      os_makefile << "WCS_INCLUDE_DIR = " + std::string(WCS_INCLUDE_DIR) + "\n";
      os_makefile << "LIBRARY_FLAGS = " + std::string(CMAKE_CXX_SHARED_LIBRARY_FLAGS) + "\n\n";
      os_makefile << "all: " + shared_lib + "\n\n";
    }

    // commands to build object file for each source file
    for (size_t i = 2u; i < m_ostreams.size(); ++i)
    {
      const std::string& src_filename
        = get_subpath(m_tmp_dir, m_ostreams[i].first);

      // This block updates no state of the current object. It only generates
      // file I/O.
      if (src_filename.empty()) {
        WCS_THROW("\n No source file to compile! Run generate_code() first.");
        continue;
      }

      std::string dir, stem, ext;

      extract_file_component(src_filename, dir, stem, ext);

      const std::string obj_filename = stem + ".o";
      obj_files += ' ' + obj_filename;

      const std::string suppress_warnings = " -Wno-unused ";
      const std::string compilation_log = (m_save_log? " 2>> wcs_jit_log.txt" : "");
      const std::string tmp_file = (m_cleanup? src_filename : "");

      std::string cmd1
        = std::string("$(CXX) $(CXXFLAGS)") + suppress_warnings
        + " -fPIC $(WCS_INCLUDE_DIR) "
        + " -c " + src_filename + compilation_log;

      os_makefile << obj_filename + ": " + src_filename + ' ' + hdr_filename + "\n"
                   + "\t" + cmd1 + "\n\n";

      //int ret = build(cmd1, obj_filename, tmp_file, compilation_log);
    }

    { // command for final linking
      std::string cmd2
        = std::string("$(CXX) $(CXXFLAGS)")
        + " -fPIC $(LIBRARY_FLAGS) "
        + " -shared -Wl,--export-dynamic " + obj_files + " -o " + m_lib_filename;

      os_makefile << shared_lib + ": " + obj_files + "\n"
                   + "\t" + cmd2 + "\n\n";
      os_makefile << "clean: \n\t@rm -f " + obj_files + " " + m_lib_filename + "\n";

      //int ret = build(cmd2, m_lib_filename, obj_files, "");
    }
  }

  close_ostream(m_ostreams[0].second);
  m_regen = false; // finished generating Makefile
 #if defined(_OPENMP)
  #pragma omp master
 #endif // defined(_OPENMP)
  {
    sync_directory(m_tmp_dir);
  }
 #if defined(_OPENMP)
  #pragma omp barrier
 #endif // defined(_OPENMP)

#if 0
  std::ifstream makefile(m_ostreams[0].first);
  if (makefile.is_open()) {
    std::cout << "Reading " + m_ostreams[0].first << std::endl;
    std::cout << makefile.rdbuf() << std::endl;
  } else {
    std::cout << "Cannot read " + m_ostreams[0].first << std::endl;
  }
#endif

  return obj_files;
}

std::string generate_cxx_code::compile_code()
{
 #if defined(_OPENMP)
  #pragma omp master
 #endif // defined(_OPENMP)
  {
    std::cerr << "JIT compiling ..." << std::endl;
  }
  int ret = EXIT_SUCCESS;

  if (!m_regen) {
    close_ostream(m_ostreams[0].second);
    return m_lib_filename;
  }

  if (m_ostreams.size() < 4u || m_ostreams[0].first.empty() ||
      m_ostreams[1].first.empty() || m_ostreams[2].first.empty()) {
    WCS_THROW("\n No source file to compile! Run generate_code() first.");
    return "";
  }

  std::string obj_files = gen_makefile();

 #if defined(_OPENMP)
  #pragma omp master
 #endif // defined(_OPENMP)
  {
    unsigned parallel_compile
      = static_cast<unsigned>(std::thread::hardware_concurrency()*0.5);
    parallel_compile = std::max(1u, parallel_compile);
    parallel_compile
      = std::min(parallel_compile, static_cast<unsigned>(m_ostreams.size()-2));

    if (m_num_compiling_threads != 0u) {
      parallel_compile
        = std::min(parallel_compile, m_num_compiling_threads);
    }

    std::string cmd3 = "pushd " + m_tmp_dir
                     + "; make -j " + std::to_string(parallel_compile)
                     + " -f " + m_ostreams[0].first + " all; popd";
    std::cout << cmd3 << std::endl;
    ret = build(cmd3, m_lib_filename, obj_files, "");
  }

  if (ret == EXIT_FAILURE) return "";

  return m_lib_filename;
}

std::vector<std::string> generate_cxx_code::get_src_filenames() const
{
  std::vector<std::string> src_filenames;
  src_filenames.reserve(m_ostreams.size());
  for (size_t i = 0u; i < m_ostreams.size(); ++i) {
    src_filenames.emplace_back(m_ostreams[i].first);
  }
  return src_filenames;
}

std::string generate_cxx_code::get_lib_filename() const
{
  return m_lib_filename;
}

} // end of namespace wcs

#endif // defined(WCS_HAS_SBML)
