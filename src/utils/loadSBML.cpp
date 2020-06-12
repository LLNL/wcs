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
  /*Konstantia*/
  Model* modelk= document->getModel();
  ListOfCompartments* compartmentsk= modelk->getListOfCompartments();
  int compSize=compartmentsk->size();
  ListOfSpecies* speciesk=modelk->getListOfSpecies();
  int speciesSize=speciesk->size();
  ListOfReactions* reactionsk=modelk->getListOfReactions();
  int reactionsSize=reactionsk->size();
  ListOfParameters* parametersk=modelk->getListOfParameters();
  int parametersSize=parametersk->size();  
  ListOfUnitDefinitions* unitsk=modelk->getListOfUnitDefinitions();
  int unitsSize=unitsk->size();
  /*Konstantia*/


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
	      /*Konstantia*/	   
	      << "\n\tNum species= "<< modelk->getNumSpecies() /*Konstantia*/
	      << "\tNum reactions= "<< modelk->getNumReactions() /*Konstantia*/
	      
	      << "\n\n//Created by WCS"
	      << "\nmodel *"<< modelk->getId()<<"()"  /*getName or getID*/
	      << "\n\n  // Compartments and Species:";  //<<compSize;
	      for(int ic=0; ic<compSize; ic++){
	      	std::cout << "\n  compartment "<< compartmentsk->get(ic)->getIdAttribute()<<";\n";  /*getName or getID*/
	      }
	      std::cout<<"  species";
	      int cnt=0;
	      for(int ic=0; ic<speciesSize; ic++){
		if (!speciesk->get(ic)->getHasOnlySubstanceUnits()){      
                  if (cnt==0){	            		
                     std::cout <<" "<<  speciesk->get(ic)->getIdAttribute()<<" in "<<speciesk->get(ic)->getCompartment();  
		  }else{	             	
                     std::cout <<", "<<  speciesk->get(ic)->getIdAttribute()<<" in "<<speciesk->get(ic)->getCompartment();  
		  }
		  cnt++;
		}
	      }
	      cnt=0;
	      for(int ic=0; ic<speciesSize; ic++){
                if (speciesk->get(ic)->getHasOnlySubstanceUnits()){
                  if (cnt==0){ 
                     std::cout<<";\n  substanceOnly species";
                     std::cout <<" $"<<  speciesk->get(ic)->getIdAttribute()<<" in "<<speciesk->get(ic)->getCompartment();
                  }else{                        
                     std::cout <<", $"<<  speciesk->get(ic)->getIdAttribute()<<" in "<<speciesk->get(ic)->getCompartment();
                  }
		  cnt++;
                }
              }
	      std::cout<<";\n\n  // Reactions:";
	      for(int ic=0; ic<reactionsSize; ic++){
                std::cout <<"\n  "<<  reactionsk->get(ic)->getIdAttribute()<<": ";
		int reactSize=reactionsk->get(ic)->getNumReactants();
		int prodSize=reactionsk->get(ic)->getNumProducts();
		//std::cout<<"\n "<<reactSize<<" "<<prodSize;
		for(int ire=0;ire<reactSize;ire++)
	        {
	          if (ire!=0){ 
                    std::cout<<" + ";
		  }
		  if (speciesk->get(reactionsk->get(ic)->getReactant(ire)->getSpecies())->getHasOnlySubstanceUnits()){
		    std::cout<<"$";
		  }			  
	          std::cout << reactionsk->get(ic)->getReactant(ire)->getSpecies();
		}
	      	std::cout<<" => "; 
		for(int ire=0;ire<prodSize;ire++)
		{              
       		  if (ire!=0){
                    std::cout<<" + ";
                  }
	          if (speciesk->get(reactionsk->get(ic)->getProduct(ire)->getSpecies())->getHasOnlySubstanceUnits()){
                    std::cout<<"$";
                  }

	          std::cout << reactionsk->get(ic)->getProduct(ire)->getSpecies()<<"; ";
		}
		std::cout << reactionsk->get(ic)->getKineticLaw()->getFormula()<<"; ";
              }
              std::cout<< "\n\n  // Species initializations:";  
              //cnt=0;
              for(int ic=0; ic<speciesSize; ic++){
		  if(speciesk->get(ic)->getInitialAmount()==0){
			   std::cout <<"\n  "<<  speciesk->get(ic)->getIdAttribute()<<" = 0";
		  }else{
		    if(!isnan(speciesk->get(ic)->getInitialAmount())){
                     std::cout <<"\n  "<< speciesk->get(ic)->getIdAttribute()<<" = "<<speciesk->get(ic)->getInitialAmount()<<"/"<<speciesk->get(ic)->getCompartment()<<";";
		    }else if(!isnan(speciesk->get(ic)->getInitialConcentration())){
			    if (speciesk->get(ic)->getInitialConcentration()==0){
			       std::cout <<"\n  "<< speciesk->get(ic)->getIdAttribute()<<" = 0;";
			    }else{
                               std::cout <<"\n  "<< speciesk->get(ic)->getIdAttribute()<<" = "<<speciesk->get(ic)->getInitialConcentration()<<"/"<<speciesk->get(ic)->getCompartment()<<";";
			    } 
		    }		    
                  }
              }
	      std::cout<< "\n\n  // Compartments initializations:";  
              for(int ic=0; ic<compSize; ic++){
                std::cout << "\n  "<< compartmentsk->get(ic)->getIdAttribute()<<" = "<<compartmentsk->get(ic)->getSize() <<";"; 
                if (compartmentsk->get(ic)->getUnits()!=""){
		   std::cout << "\n  "<< compartmentsk->get(ic)->getIdAttribute()<<" has "<<compartmentsk->get(ic)->getUnits() <<";"; 
		}
	      }

               std::cout<< "\n\n  // Variable initializations:";  
              for(int ic=0; ic<parametersSize; ic++){
                std::cout << "\n  "<< parametersk->get(ic)->getIdAttribute()<<" = "<<parametersk->get(ic)->getValue() <<";";       
                if (parametersk->get(ic)->getUnits()!=""){
		   std::cout << "\n  "<< parametersk->get(ic)->getIdAttribute()<<" has "<<parametersk->get(ic)->getUnits() <<";";         
                }
	      }
               
              std::cout<< "\n\n  // Other declarations:"<<"\n  const ";
              for(int ic=0; ic<compSize; ic++){
                std::cout << compartmentsk->get(ic)->getIdAttribute();
              }
              for(int ic=0; ic<parametersSize; ic++){
                std::cout <<", "<< parametersk->get(ic)->getIdAttribute();   
              }
	      std::cout <<";";

              std::cout<< "\n\n  // Unit definitions:";
              for(int ic=0; ic<unitsSize; ic++){
                std::cout << "\n  unit "<< unitsk->get(ic)->getIdAttribute();
		if (unitsk->get(ic)->getIdAttribute()=="time"){
			std::cout<<"_unit";
		}
		std::cout<<" = ";
		ListOfUnits* unitslistk=unitsk->get(ic)->getListOfUnits();
		int unitslistSize=unitslistk->size();
		for(int iu=0; iu<unitslistSize;iu++){
                  std::cout<<UnitKind_toString(unitslistk->get(iu)->getKind());	
                  if (unitslistk->get(iu)->getExponent()>1){
			  std::cout<<"^"<<unitslistk->get(iu)->getExponent();
		  }
		  std::cout<<";";

		}
              }

	      std::cout<< "\n\n  // Display names:";
              for(int ic=0; ic<unitsSize; ic++){
                if (unitsk->get(ic)->getIdAttribute()!=unitsk->get(ic)->getName()&&unitsk->get(ic)->getName()!=""){
                   std::cout << "\n  "<< unitsk->get(ic)->getIdAttribute()<<" is \"" <<unitsk->get(ic)->getName() <<"\";";
                }

	        if (unitsk->get(ic)->getIdAttribute()=="time"){
                 // if(){        
			std::cout<<"\n  time_unit is \""<<unitsk->get(ic)->getName()<<"\";";
		//}
		}

              }

              for(int ic=0; ic<parametersSize; ic++){
                if (parametersk->get(ic)->getIdAttribute()!=parametersk->get(ic)->getName()&&parametersk->get(ic)->getName()!=""){
                   std::cout << "\n  "<< parametersk->get(ic)->getIdAttribute()<<" is \"" <<parametersk->get(ic)->getName() <<"\";";
                }
              }

	    std::cout<<"\nend\n";
 
	     
	      /*Konstantia*/
              
	      //<< std::endl;
  }

  Model* model = document->getModel();

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
