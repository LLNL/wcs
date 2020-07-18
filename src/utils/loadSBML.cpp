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

#include "sbml_utils.hpp"

#if defined(WCS_HAS_SBML)
#include <iostream>

#include <sbml/SBMLTypes.h>
#include <sbml/common/extern.h>


LIBSBML_CPP_NAMESPACE_USE

int main(int argc, char** argv)
{
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " filename" << std::endl;
    return 0;
  }

  const char* filename = argv[1];
  SBMLReader reader;

  SBMLDocument* document = reader.readSBML(filename);
  if (document == nullptr) {
    std::cout << "Faile to get SBML document from " << filename << std::endl;
    return -1;
  }

  Model* model = document->getModel();
  if (model == nullptr) {
    std::cout << "Faile to get SBML model from " << filename << std::endl;
    delete document;
    return -1;
  }

  ListOfCompartments* compartment_list = model->getListOfCompartments();
  const unsigned int num_compartments = compartment_list->size();
  ListOfSpecies* species_list = model->getListOfSpecies();
  const unsigned int num_species = species_list->size();
  ListOfReactions* reaction_list = model->getListOfReactions();
  const unsigned int num_reactions = reaction_list->size();
  ListOfParameters* parameter_list = model->getListOfParameters();
  const unsigned int num_parameters = parameter_list->size();
  ListOfUnitDefinitions* unit_list = model->getListOfUnitDefinitions();
  const unsigned int num_units = unit_list->size();


  const unsigned num_errors = document->getNumErrors();

  if (num_errors > 0u) {
    std::cout << num_errors << " error(s) in reading "
              << filename << std::endl;

    document->printErrors(std::cerr);

    delete document;
    return static_cast<int>(num_errors);
  }

  std::cout << "SBML file " << filename << ":" << std::endl
            << "\tLevel = " << document->getLevel()
            << "\tVersion = " << document->getVersion()
            << "\n\tNum species = " << model->getNumSpecies()
            << "\tNum reactions = " << model->getNumReactions()
            << "\n\n//Created by WCS"
            << "\nmodel *" << model->getId() << "()"  ///getName or getID

            // Print out Compartments and Species
            << "\n\n  // Compartments and Species:";
  for (unsigned int ic = 0u; ic<num_compartments; ic++) {
    std::cout << "\n  compartment " << compartment_list->get(ic)->getIdAttribute()
    << ";\n";  ///getName or getID
  }
  std::cout << "  species";
  unsigned int cnt = 0u;
  for (unsigned int ic = 0u; ic < num_species; ic++) {
    const LIBSBML_CPP_NAMESPACE::Species& species = *(species_list->get(ic));
    if (!species.getHasOnlySubstanceUnits()){
      if (cnt == 0u){
        std::cout << " " <<  species.getIdAttribute()
                  << " in "<< species.getCompartment();
             } else {
        std::cout << ", " <<  species.getIdAttribute()
                  << " in "<< species.getCompartment();
      }
      cnt++;
    }
  }
  cnt = 0u;
  for (unsigned int ic = 0u; ic < num_species; ic++) {
    const LIBSBML_CPP_NAMESPACE::Species& species = *(species_list->get(ic));
    if (species.getHasOnlySubstanceUnits()){
      if (cnt == 0u) {
        std::cout << ";\n  substanceOnly species";
        std::cout << " $" <<  species.getIdAttribute()
                  << " in "<< species.getCompartment();
      } else {
        std::cout << ", $" << species.getIdAttribute()
                  << " in "<< species.getCompartment();
      }
      cnt++;
    }
  }

  // Print out Reactions
  std::cout << ";\n\n  // Reactions:";
  for (unsigned int ic = 0u; ic < num_reactions; ic++) {
    const LIBSBML_CPP_NAMESPACE::Reaction& reaction = *(reaction_list->get(ic));
    std::cout << "\n  " <<  reaction.getIdAttribute() << ": ";
    unsigned int num_reactants = reaction.getNumReactants();
    unsigned int num_products = reaction.getNumProducts();
    for (unsigned int ire = 0u; ire < num_reactants; ire++){
      if (ire!=0) {
        std::cout << " + ";
      }
      if (species_list->get(reaction.getReactant(ire)->getSpecies())->
          getHasOnlySubstanceUnits())
      {
        std::cout << "$";
      }
      std::cout << reaction.getReactant(ire)->getSpecies();
    }
    std::cout << " => ";
    for (unsigned int ire = 0u; ire < num_products; ire++) {
      if (ire!=0) {
        std::cout << " + ";
      }
      if (species_list->get(reaction.getProduct(ire)->getSpecies())->
          getHasOnlySubstanceUnits())
      {
        std::cout << "$";
      }
      std::cout << reaction.getProduct(ire)->getSpecies() << "; ";
    }
    std::cout << reaction.getKineticLaw()->getFormula() << "; ";
    //std::cout << "\n "<< SBML_formulaToString(reaction->getKineticLaw()->getMath())<<"; ";

    std::string formula
      = SBML_formulaToString(reaction.getKineticLaw()->getMath()); //char *
    const std::string toReplace("pow(");
    const std::string toReplace2(", ");
    const std::string toReplace3(")");
    std::string wholeformula("");
    size_t pos = formula.find(toReplace);

    while (pos < formula.length()) {
      if (pos != std::string::npos) {
        formula.replace(pos,toReplace.length(),"");
      }
      size_t pos2 = formula.find(toReplace2,pos);
      formula.replace(pos2,toReplace2.length(),"^");
      size_t pos3 = formula.find(toReplace3,pos);
      formula.replace(pos3,toReplace3.length(),"");
      pos = formula.find(toReplace);
    }
    std::cout << "\n  ";

    for (unsigned int ic = 0u; ic < num_parameters; ic++) {
      const LIBSBML_CPP_NAMESPACE::Parameter& parameter = *(parameter_list->get(ic));
      std::string toFindPar(parameter.getIdAttribute());
      size_t posPar = formula.find(toFindPar);
      std::string parametervalue = std::to_string(parameter.getValue()).substr(0,
        std::to_string(parameter.getValue()).find(".") + 3);

      if (posPar != std::string::npos) {
        wholeformula = wholeformula + "var " + parameter.getIdAttribute()
                     + " := " + parametervalue +  "; ";
        std::cout << "var " << parameter.getIdAttribute()
                  << " := " << parameter.getValue() << "; ";
      }
    }
    std::cout << "m_rate := " << formula << ";\n";
    wholeformula = wholeformula + "m_rate := " + formula + ";\n";
    //std::cout << wholeformula;

    //const ASTNode*  math = reaction_list->get(ic)->getKineticLaw()->getMath();
    //wcs::sbml_utils sbml_o;
    using reaction_parameters = std::unordered_set<std::string>;
    reaction_parameters pset;

    //pset= sbml_o.get_symbol_table_of_formula(*math);
    //pset = wcs::sbml_utils::get_symbol_table_of_formula(*math);
  }

  // Print out Species initializations
  std::cout << "\n\n  // Species initializations:";
  for (unsigned int ic = 0u; ic < num_species; ic++) {
    const LIBSBML_CPP_NAMESPACE::Species& species = *(species_list->get(ic));
    if (species.getInitialAmount() == 0) {
      std::cout << "\n  " << species.getIdAttribute() << " = 0";
    } else {
      if (!isnan(species.getInitialAmount())) {
        std::cout << "\n  " << species.getIdAttribute()
                  << " = " << species.getInitialAmount() << ";";
      } else if (!isnan(species.getInitialConcentration())) {
        if (species.getInitialConcentration() == 0) {
          std::cout << "\n  " << species.getIdAttribute() << " = 0;";
        } else {
          std::cout << "\n  " << species.getIdAttribute()
                    << " = " << species.getInitialConcentration() << ";";
        }
      }
    }
  }

  // Print out Compartments initializations
  std::cout << "\n\n  // Compartments initializations:";
  for (unsigned int ic = 0u; ic < num_compartments; ic++) {
    const LIBSBML_CPP_NAMESPACE::Compartment& compartment
      = *(compartment_list->get(ic));
    std::cout << "\n  " << compartment.getIdAttribute()
              << " = " << compartment.getSize() << ";";
    if (compartment.getUnits() != "") {
      std::cout << "\n  " << compartment.getIdAttribute()
                << " has " << compartment.getUnits() << ";";
    }
  }

  // Print out Variable initializations
  std::cout << "\n\n  // Variable initializations:";
  for (unsigned int ic = 0u; ic < num_parameters; ic++) {
    const LIBSBML_CPP_NAMESPACE::Parameter& parameter
      = *(parameter_list->get(ic));
    std::cout << "\n  " << parameter.getIdAttribute()
              << " = " << parameter.getValue() << ";";
    if (parameter.getUnits() != "") {
      std::cout << "\n  " << parameter.getIdAttribute()
                << " has " << parameter.getUnits() << ";";
    }
  }

  // Print out other declarations
  std::cout << "\n\n  // Other declarations:" << "\n  const ";
  for (unsigned int ic = 0u; ic < num_compartments; ic++) {
    std::cout << compartment_list->get(ic)->getIdAttribute();
  }
  for (unsigned int ic = 0u; ic < num_parameters; ic++) {
    std::cout << ", " << parameter_list->get(ic)->getIdAttribute();
  }
  std::cout << ";";

  // Print out Unit definitions
  std::cout << "\n\n  // Unit definitions:";
  for (unsigned int ic = 0u; ic < num_units; ic++) {
    const LIBSBML_CPP_NAMESPACE::UnitDefinition& unit = *(unit_list->get(ic));
    std::cout << "\n  unit "<< unit.getIdAttribute();
    if (unit.getIdAttribute() == "time") {
      std::cout << "_unit";
    }
    std::cout << " = ";
    const ListOfUnits* subunit_list = unit.getListOfUnits();
    const unsigned int subnum_units = subunit_list->size();
    //std::cout<<"\n "<<unitslistSize<<"\n";
    for (unsigned int iu = 0u; iu < subnum_units; iu++) {
      const LIBSBML_CPP_NAMESPACE::Unit& subunit = *(subunit_list->get(iu));
      if (subunit.getExponent() == -1) {
        std::cout << "/";
      }
      std::cout << UnitKind_toString(subunit.getKind());
      if (subunit.getExponent() > 1) {
        std::cout << "^" << subunit.getExponent();
      }
      if (iu == subnum_units-1) {
        std::cout << ";";
      }
    }
  }

  // Print out Display names
  std::cout << "\n\n  // Display names:";
  for (unsigned int ic = 0u; ic < num_units; ic++) {
    const LIBSBML_CPP_NAMESPACE::UnitDefinition& unit = *(unit_list->get(ic));
    if (unit.getIdAttribute() != unit.getName() && unit.getName() != "") {
      std::cout << "\n  " << unit.getIdAttribute()
                << " is \"" << unit.getName() << "\";";
    }

    if (unit.getIdAttribute() == "time") {
      std::cout << "\n  time_unit is \"" << unit.getName() << "\";";
    }
  }

  for (unsigned int ic = 0u; ic < num_parameters; ic++) {
    const LIBSBML_CPP_NAMESPACE::Parameter& parameter
      = *(parameter_list->get(ic));
    if (parameter.getIdAttribute() != parameter.getName() &&
        parameter.getName() != "")
    {
      std::cout << "\n  " << parameter.getIdAttribute()
                << " is \"" << parameter.getName() << "\";";
    }
  }

  std::cout << "\nend\n";

  delete document;
  return 0;
}

#endif // defined(WCS_HAS_SBML)
