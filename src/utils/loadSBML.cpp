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

#if defined(WCS_HAS_SBML)
#include <iostream>

#include <sbml/SBMLTypes.h>
#include <sbml/common/extern.h>


LIBSBML_CPP_NAMESPACE_USE

BEGIN_C_DECLS

int main(int argc, char** argv)
{
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " filename" << std::endl;
    return 0;
  }

  const char* filename = argv[1];
  SBMLReader reader;
  SBMLDocument* document = reader.readSBML(filename);
  Model* model = document->getModel();
  ListOfCompartments* compartmentslist = model->getListOfCompartments();
  int compartmentsSize = compartmentslist->size();
  ListOfSpecies* specieslist = model->getListOfSpecies();
  int speciesSize = specieslist->size();
  ListOfReactions* reactionslist = model->getListOfReactions();
  int reactionsSize = reactionslist->size();
  ListOfParameters* parameterslist = model->getListOfParameters();
  int parametersSize = parameterslist->size();  
  ListOfUnitDefinitions* unitslist = model->getListOfUnitDefinitions();
  int unitsSize = unitslist->size();


  const unsigned num_errors = document->getNumErrors();

  if (num_errors > 0u) {
    std::cout << num_errors << " error(s) in reading "
              << filename << std::endl;

    document->printErrors(std::cerr);

    delete document;
    return static_cast<int>(num_errors);
  } else {
    std::cout << "SBML file " << filename << ":" << std::endl
              << "\tLevel = " << document->getLevel()
              << "\tVersion = " << document->getVersion()
	      << "\n\tNum species = " << model->getNumSpecies() 
	      << "\tNum reactions = " << model->getNumReactions() 
	      
	      << "\n\n//Created by WCS"
	      << "\nmodel *" << model->getId() << "()"  ///getName or getID

	      /// Print out Compartments and Species
	      << "\n\n  // Compartments and Species:";  
    for (int ic = 0; ic<compartmentsSize; ic++) {
      std::cout << "\n  compartment " << compartmentslist->get(ic)->getIdAttribute() << ";\n";  ///getName or getID
    } 
    std::cout << "  species";
    int cnt = 0;
    for (int ic = 0; ic < speciesSize; ic++) {
      if (!specieslist->get(ic)->getHasOnlySubstanceUnits()){      
        if (cnt == 0){	            		
          std::cout << " " <<  specieslist->get(ic)->getIdAttribute() << " in "<<specieslist->get(ic)->getCompartment();  
       	} else {	             	
          std::cout << ", " <<  specieslist->get(ic)->getIdAttribute() << " in "<<specieslist->get(ic)->getCompartment();  
        }  
        cnt++;
      }
    }
    cnt = 0;
    for (int ic = 0; ic < speciesSize; ic++) {
      if (specieslist->get(ic)->getHasOnlySubstanceUnits()){
        if (cnt == 0) { 
          std::cout << ";\n  substanceOnly species";
          std::cout << " $" <<  specieslist->get(ic)->getIdAttribute() << " in "<<specieslist->get(ic)->getCompartment();
        } else {                        
          std::cout << ", $" <<  specieslist->get(ic)->getIdAttribute() << " in "<<specieslist->get(ic)->getCompartment();
	}	
      	cnt++;
      }
    }
  
    /// Print out Reactions
    std::cout << ";\n\n  // Reactions:";
    for (int ic = 0; ic < reactionsSize; ic++) {
      std::cout << "\n  " <<  reactionslist->get(ic)->getIdAttribute() << ": ";
      int reactSize = reactionslist->get(ic)->getNumReactants();
      int prodSize = reactionslist->get(ic)->getNumProducts();
      for (int ire = 0; ire < reactSize; ire++){
        if (ire!=0) {  
          std::cout << " + ";
        }
	if (specieslist->get(reactionslist->get(ic)->getReactant(ire)->getSpecies())->getHasOnlySubstanceUnits()){
          std::cout << "$";
	}			  
        std::cout << reactionslist->get(ic)->getReactant(ire)->getSpecies();
      }
      std::cout << " => "; 
      for (int ire = 0; ire < prodSize; ire++) {              
        if (ire!=0) { 
          std::cout << " + ";
        }
        if (specieslist->get(reactionslist->get(ic)->getProduct(ire)->getSpecies())->getHasOnlySubstanceUnits()) {
          std::cout << "$";
        }
        std::cout << reactionslist->get(ic)->getProduct(ire)->getSpecies() << "; ";
      }
      std::cout << reactionslist->get(ic)->getKineticLaw()->getFormula() << "; ";
      ///std::cout << "\n "<< SBML_formulaToString( reactionslist->get(ic)->getKineticLaw()->getMath())<<"; ";
               
      std::string formula = SBML_formulaToString(reactionslist->get(ic)->getKineticLaw()->getMath()); //char *
      std::string toReplace("pow(");
      std::string toReplace2(", ");
      std::string toReplace3(")");
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
      for (int ic = 0; ic < parametersSize; ic++) {
        std::string toFindPar(parameterslist->get(ic)->getIdAttribute());
        size_t posPar = formula.find(toFindPar);
	std::string parametervalue = std::to_string(parameterslist->get(ic)->getValue()).substr(0, std::to_string(parameterslist->get(ic)->getValue()).find(".") + 3);
        if (posPar != std::string::npos) {
	  wholeformula = wholeformula + "var " + parameterslist->get(ic)->getIdAttribute() + " := " + parametervalue +  "; ";	
          std::cout << "var " <<parameterslist->get(ic)->getIdAttribute() << " := " << parameterslist->get(ic)->getValue() << "; ";
        } 
      }
      std::cout << "m_rate := " << formula << ";\n";
      wholeformula = wholeformula + "m_rate := " + formula + ";\n";
      ///std::cout << wholeformula;

      /**std::cout << "FunctionDefinition " << reactionsk->get(ic)->getIdAttribute();
      const ASTNode*  math = reactionsk->get(ic)->getKineticLaw()->getMath();   
      if (math->getNumChildren() > 1) {
        std::cout<<" " << math->getNumChildren()<<"\n";
	std::cout << "\n(" << ( math->getLeftChild() )->getName();
        for (int n = 0; n < math->getNumChildren() - 1; ++n){
       	  std::cout <<", " << ( math->getChild(n) )->getName();
        }
      }   

      std::cout <<") := ";
      if (math->getNumChildren() == 0){ 
        std::cout << "(no body defined)";
      } else {
        math    = math->getChild(math->getNumChildren() - 1);
        char* formula = SBML_formulaToString(math);
        std::cout << formula ;
      }*/

   
    }  
   
    /// Print out Species initializations
    std::cout << "\n\n  // Species initializations:";  
    for (int ic = 0; ic < speciesSize; ic++) {
      if (specieslist->get(ic)->getInitialAmount() == 0) {
        std::cout << "\n  " <<  specieslist->get(ic)->getIdAttribute() << " = 0";
      } else {
        if (!isnan(specieslist->get(ic)->getInitialAmount())){
          std::cout << "\n  " << specieslist->get(ic)->getIdAttribute() << " = " << specieslist->get(ic)->getInitialAmount() << ";"; /// << "/" <<specieslist->get(ic)->getCompartment() << ";";
       	} else if (!isnan(specieslist->get(ic)->getInitialConcentration())) {
          if (specieslist->get(ic)->getInitialConcentration() == 0) { 
	    std::cout << "\n  " << specieslist->get(ic)->getIdAttribute() << " = 0;";
	  } else {
            std::cout << "\n  " << specieslist->get(ic)->getIdAttribute() << " = " << specieslist->get(ic)->getInitialConcentration() << ";"; /// << "/" << specieslist->get(ic)->getCompartment() << ";";
	  } 
	}		
      }
    }
    
    /// Print out Compartments initializations
    std::cout << "\n\n  // Compartments initializations:";  
    for (int ic = 0; ic < compartmentsSize; ic++) {
      std::cout << "\n  " << compartmentslist->get(ic)->getIdAttribute() << " = " <<compartmentslist->get(ic)->getSize() << ";"; 
      if (compartmentslist->get(ic)->getUnits() != "") {
        std::cout << "\n  " << compartmentslist->get(ic)->getIdAttribute() << " has " << compartmentslist->get(ic)->getUnits() << ";"; 
      }
    }

    /// Print out Variable initializations
    std::cout << "\n\n  // Variable initializations:";  
    for (int ic = 0; ic < parametersSize; ic++) {
      std::cout << "\n  " << parameterslist->get(ic)->getIdAttribute() << " = " << parameterslist->get(ic)->getValue() << ";";       
      if (parameterslist->get(ic)->getUnits() != "") {
        std::cout << "\n  " << parameterslist->get(ic)->getIdAttribute() << " has " << parameterslist->get(ic)->getUnits() << ";";         
      }
    }
           
    /// Print out other declarations
    std::cout << "\n\n  // Other declarations:" << "\n  const ";
    for (int ic = 0; ic < compartmentsSize; ic++) {
      std::cout << compartmentslist->get(ic)->getIdAttribute();
    }
    for (int ic = 0; ic < parametersSize; ic++) {
      std::cout << ", " << parameterslist->get(ic)->getIdAttribute();   
    }
    std::cout << ";";

    /// Print out Unit definitions
    std::cout << "\n\n  // Unit definitions:";
    for (int ic = 0; ic < unitsSize; ic++) {
      std::cout << "\n  unit "<< unitslist->get(ic)->getIdAttribute();
      if (unitslist->get(ic)->getIdAttribute() == "time") {
        std::cout << "_unit";
      }
      std::cout << " = ";
      ListOfUnits* subunitslist = unitslist->get(ic)->getListOfUnits();
      int subunitsSize = subunitslist->size();
      //std::cout<<"\n "<<unitslistSize<<"\n";
      for (int iu = 0; iu < subunitsSize ; iu++) {
        if (subunitslist->get(iu)->getExponent() == -1) {
          std::cout << "/";
        }
	std::cout << UnitKind_toString(subunitslist->get(iu)->getKind());	
        if (subunitslist->get(iu)->getExponent() > 1) {
          std::cout << "^" << subunitslist->get(iu)->getExponent();
	}
	if (iu == subunitsSize-1) {
	  std::cout << ";";
	}
      }
    }

    /// Print out Display names
    std::cout << "\n\n  // Display names:";
    for (int ic = 0; ic < unitsSize; ic++) {
      if (unitslist->get(ic)->getIdAttribute() != unitslist->get(ic)->getName() && unitslist->get(ic)->getName() != "") {
        std::cout << "\n  " << unitslist->get(ic)->getIdAttribute() << " is \"" << unitslist->get(ic)->getName() << "\";";
      }

      if (unitslist->get(ic)->getIdAttribute() == "time") {
        std::cout << "\n  time_unit is \"" << unitslist->get(ic)->getName() << "\";";
      }
    }

    for (int ic = 0; ic < parametersSize; ic++) {
      if (parameterslist->get(ic)->getIdAttribute() != parameterslist->get(ic)->getName() && parameterslist->get(ic)->getName() != ""){
        std::cout << "\n  " << parameterslist->get(ic)->getIdAttribute() << " is \"" << parameterslist->get(ic)->getName() << "\";";
      }
    }

    std::cout << "\nend\n";
  }


  if (model == NULL) {
    std::cout << "Faile to get model from " << filename << std::endl;
    delete document;
    return -1;
  }

  if (model->isSetSBOTerm()) {
    std::cout << "Num species: " << model->getNumSpecies() << std::endl;
    std::cout << "Num reactions: " << model->getNumReactions() << std::endl;
  }

  delete document;
  return 0;
}

END_C_DECLS
#endif // defined(WCS_HAS_SBML)
