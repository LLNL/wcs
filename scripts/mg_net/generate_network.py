import sys, argparse
import os.path
#import random
import numpy
from argparse import RawTextHelpFormatter
from libsbml import *

import xml.etree.ElementTree as ET

#import xml.dom.minidom as md

def main (args):
  """Usage: printNotes filename
  """

  parser = argparse.ArgumentParser(description='''Generate subnetworks based on the network of the Mycoplasma genitalium organism (mg_net network).
See README.md file for more information.''', formatter_class=RawTextHelpFormatter)
  parser._action_groups.pop()
  requiredNamed = parser.add_argument_group('required arguments and inputs')
  optionalNamed = parser.add_argument_group('optional arguments')
  requiredNamed.add_argument('-s', metavar='<sbml_filename>', nargs=1, help='Specify the SBML filename of mg_net network', required=True, dest="sbml_filename")
  requiredNamed.add_argument('-g', metavar='<graphML_filename>', nargs=1, help='Specify the GraphML filename of mg_net network', required=True, dest="graphml_filename")
  requiredNamed.add_argument('-i', metavar='<input_count_species_filename>', nargs=1, help='''Specify the SBML filename for the initial count of the species of the metabolism cellular 
process contained in the specific network''', dest="species_values_filename")
  optionalNamed.add_argument('-o', metavar='<output_filename>', nargs=1, help='Specify the output filename', dest="output_filename")
  optionalNamed.add_argument('-x', metavar='<approach_id>', type=int, nargs=1, help='''Specify which approach you want to generate (1 for the first approach or 2 for the second approach
or 3 for the third approach or 23 for the combined approach of 2 and 3 )''', dest="approach_x")
  optionalNamed.add_argument('-p', metavar='<process_id>', nargs='+', help= '''Specify the cellular processes to be included in the new subnetwork (default: 
all the processes) by listing the process_ids separated by a space.:

1. Replication
2. Transcription
3. Protein-DNA Interaction_
4. RNA Degradation
5. Translation
6. Protein Degradation
7. Metabolism
8. Protein Maturation
9. Macromolecular Complexation
10. Protein-DNA Interaction
11. tRNA Aminoacylation
12. Terminal Organelle Assembly
13. RNA Processing
14. Replication Initiation
15. Cellular Division
  ''', dest="proc", default=list(range(1, 15+1)))
  optionalNamed.add_argument('-u', metavar=('<min>', '<max>'), type=int, nargs=2, help='''Specify the range of the the “discrete uniform” distribution
used in the random initialization of the initial copy-number of each 
species  (default from 30000 to 500000)''',
  dest="range_values", default=["30000", "500000"])
  optionalNamed.add_argument('-n', metavar=('<min>', '<max>'), nargs=2, type=float, help='''Specify the range of the normal distribution 
for the initialization of reaction rate coefficient k used in the
reaction rate formulas (default from 2.5 to 7.5). k is the same for all the reactions in a cellular process.''',
  dest="normal_distr", default=["2.5", "7.5"])
  optionalNamed.add_argument('-m', metavar=('<max_sum_stoichiometry>'), nargs=1, type=int, help='''Specify the maximum sum of stoichiometries as the threshold between the 
  first and third approach. Value should be between 3 and 10 (default value 7).''',
  dest="max_sum_stoich", default=["7"])


  original_stdout = sys.stdout # Save a reference to the original standard output

  try:
    args = parser.parse_args()

    #print('Range:', args.range_values)
    if args.range_values != None :
      if (int(args.range_values[0]) > int(args.range_values[1]) or int(args.range_values[0]) < 0 ) :
        print('Error: the first argument of -u has to be greater than 0 and smaller than the second argument')
        sys.exit()
      min_range_value = int(args.range_values[0])
      max_range_value = int(args.range_values[1])
      #print ("\n" + str(min_range_value) + " " + str(max_range_value))

    if args.normal_distr != None :
      if (float(args.normal_distr[0]) > float(args.normal_distr[1]) or float(args.normal_distr[0]) < 0 ) :
        print('Error: the first argument of -n has to be greater than 0 and smaller than the second argument')
        sys.exit()
      min_k_value = float(args.normal_distr[0])
      max_k_value = float(args.normal_distr[1])
    #print ("\n" + str(min_k_value) + " " + str(max_k_value))

    processes = set ()
    approach1 = True
    approach3 = True
    combined_case = False
    include_Metabolism = True
    input_used = False 
    max_sum_stoich = 7  #default value
    
    #print('Processes:', args.proc)
    if args.proc != None :
      for opt in args.proc:
        if (int(opt) > 15 or int(opt) < 0 ) :
          print('Error: the numbers corresponding to processes have to be between 1 to 15. Run \"python3 generate_network.py -h\" to see the number for each process.')
          sys.exit()
        processes.add(opt)

    #for opt in processes:
    #    print ("\n" + str(opt))

    if args.species_values_filename != None :
      input_filename = str(args.species_values_filename[0])
      document_input = readSBML(input_filename)
      if document_input.getNumErrors() > 0:
        print("Encountered the following SBML errors:" )
        document_input.printErrors()
        return 1

      model_input = document_input.getModel()

      if model_input is None:
        print("No model present." )
        return 1
      input_used = True
      include_Metabolism = False   

    
    if args.approach_x != None :
      #print ("\n" + str(args.approach_x[0])) 
      if (int(args.approach_x[0]) !=1 and int(args.approach_x[0])!=2 and int(args.approach_x[0])!=3 and int(args.approach_x[0])!=23) :
        print('Error: the approach id of -x has to be 1 or 2 or 3 or 23')
        sys.exit()
      processes.clear()  
      if (int(args.approach_x[0])==1):
        processes = {1,2,3,5,7,8,9,10,11,12,13,14,15}
        approach3 = False
      elif (int(args.approach_x[0])==2):
        processes = {4,6}
        include_Metabolism = False 
      elif (int(args.approach_x[0])==3):
        processes = {1,2,5,8,9,13}
        approach1 = False
        include_Metabolism = False  
      elif (int(args.approach_x[0])==23):
        processes2 = {4,6}
        processes = {1,2,5,8,9,13}
        combined_case = True
        approach1 = False
        include_Metabolism = False

    if args.max_sum_stoich != None :
      if (int(args.max_sum_stoich[0]) < 3 or int(args.max_sum_stoich[0]) > 10 ) :
        print('Error: the max_stoichiometry of -m has to be between 3 and 10.')
        sys.exit()
      max_sum_stoich = int(args.max_sum_stoich[0])  
      

    if args.output_filename != None :
      output_filename = str(args.output_filename[0])
      f = open(output_filename, "w")
      sys.stdout = f

  except IOError as msg:
    parser.error(str(msg))

  normal_distr = ''


  sbml_filename = str(args.sbml_filename[0])
  graphml_filename = str(args.graphml_filename[0])

  document = readSBML(sbml_filename)

  tree = ET.parse(graphml_filename)
  root = tree.getroot()

  if document.getNumErrors() > 0:
    print("Encountered the following SBML errors:" )
    document.printErrors()
    return 1

  model = document.getModel()

  if model is None:
    print("No model present." )
    return 1

  reactions_ids = set ()
  reactions_metab_ids =set ()
  reactions_dict = {}
  if len(processes) > 0 :
    processes_names = ["Replication", "Transcription","Protein-DNA Interaction ","RNA Degradation","Translation"
    ,"Protein Degradation","Metabolism","Protein Maturation","Macromolecular Complexation","Protein-DNA Interaction"
    ,"tRNA Aminoacylation","Terminal Organelle Assembly","RNA Processing","Replication Initiation","Cellular Division"]
  
  #Find the reactions for the selected processes
  #print('\n Type ='+ str(len(processes)))
  for elem in root:
    for subelem in elem:
      for subsubelem in subelem:
        for  proc in processes:
          if (subsubelem.text == processes_names[int(proc)-1]):

            reaction = model.getReaction(str(subelem.attrib.get('id'))) 
            # Delete H2O and H
            reactants = reaction.getListOfReactants()
            products = reaction.getListOfProducts()
            products.remove("c_H")
            products.remove("c_H2O")
            products.remove("e_H")
            products.remove("e_H2O")
            reactants.remove("c_H")
            reactants.remove("c_H2O")
            reactants.remove("e_H")
            reactants.remove("e_H2O")
            sum_stoich = 0
            for y in range(0, reaction.getNumReactants()):
              sum_stoich = sum_stoich + reaction.getReactant(y).getStoichiometry()
            
            if (approach3 == False and sum_stoich < max_sum_stoich) or (approach1 == False and sum_stoich >= max_sum_stoich) or (approach1 == True and approach3 == True):
              #print(processes_names[int(proc)-1])
              if processes_names[int(proc)-1] == "Metabolism" and include_Metabolism == True :
                reactions_metab_ids.add(str(subelem.attrib.get('id')))
              else:  
                reactions_ids.add(str(subelem.attrib.get('id')))
              
              reactions_dict.setdefault(int(proc), [])
              reactions_dict[int(proc)].append(str(subelem.attrib.get('id')))


  # Add the reactions for the second approach in the combined approach 2 and 3
  #print('\n Type ='+ str(len(processes)))
  if combined_case  == True:
    approach1 = True 
    for elem in root:
      for subelem in elem:
        for subsubelem in subelem:
          for  proc in processes2:
            if (subsubelem.text == processes_names[int(proc)-1]):

              reaction = model.getReaction(str(subelem.attrib.get('id'))) 
              # Delete H2O and H
              reactants = reaction.getListOfReactants()
              products = reaction.getListOfProducts()
              products.remove("c_H")
              products.remove("c_H2O")
              products.remove("e_H")
              products.remove("e_H2O")
              reactants.remove("c_H")
              reactants.remove("c_H2O")
              reactants.remove("e_H")
              reactants.remove("e_H2O")
              sum_stoich = 0
              
              if (approach1 == True and approach3 == True):
                # #print(processes_names[int(proc)-1])
                reactions_ids.add(str(subelem.attrib.get('id')))
                
                reactions_dict.setdefault(int(proc), [])
                reactions_dict[int(proc)].append(str(subelem.attrib.get('id'))) 
    approach1 = False     

 

  # To find only the species belong to the metabolism
  find_reactions_metabolism_ids = set ()
  if include_Metabolism == False :
    process_only_metabolism ={7}
    for elem in root:
      for subelem in elem:
        for subsubelem in subelem:
          for  proc in process_only_metabolism:
            if (subsubelem.text == processes_names[int(proc)-1]):

              reaction = model.getReaction(str(subelem.attrib.get('id'))) 
              # Delete H2O and H
              reactants = reaction.getListOfReactants()
              products = reaction.getListOfProducts()
              products.remove("c_H")
              products.remove("c_H2O")
              products.remove("e_H")
              products.remove("e_H2O")
              reactants.remove("c_H")
              reactants.remove("c_H2O")
              reactants.remove("e_H")
              reactants.remove("e_H2O")
              sum_stoich = 0
              
              find_reactions_metabolism_ids.add(str(subelem.attrib.get('id')))

  #print(len(find_reactions_metabolism_ids))

  #print(" Reactions:" + str(len(reactions_dict)))
  #print(" Reactions:" + str(len(reactions_ids)))
  #print(reactions_dict)

  # keys = reactions_dict.keys()
  # for x in keys:
  #   print ("\n key " + str(x) + " :")
    #rlist = reactions_dict.get(x)
    #for y in rlist:
    #  print (str(y) + ", ")

  

  species = set()
  reactions = set ()
  species_metab = set ()
  species_in_metab = set()
  #Add Metabolism species when metabolism is included in the processes
  for x in reactions_metab_ids:
    reaction = model.getReaction(x)
    
    for y in range(0, reaction.getNumReactants()):  
      species_metab.add(reaction.getReactant(y).getSpecies())

    for y in range(0, reaction.getNumProducts()):
      species_metab.add(reaction.getProduct(y).getSpecies())

    for y in range(0, reaction.getNumModifiers()):
      species_metab.add(reaction.getModifier(y).getSpecies())
 
  if include_Metabolism == False :
    #Find Metabolism species when Metabolism is not included in the processes
    for x in find_reactions_metabolism_ids:
      reaction = model.getReaction(x)
      
      for y in range(0, reaction.getNumReactants()):  
        species_in_metab.add(reaction.getReactant(y).getSpecies())

      for y in range(0, reaction.getNumProducts()):
        species_in_metab.add(reaction.getProduct(y).getSpecies())

      for y in range(0, reaction.getNumModifiers()):
        species_in_metab.add(reaction.getModifier(y).getSpecies())

  #Add the species of the rest cellular processes
  for x in reactions_ids:
    reaction = model.getReaction(x)
    
    for y in range(0, reaction.getNumReactants()):  
      if reaction.getReactant(y).getSpecies() not in species_metab:
        species.add(reaction.getReactant(y).getSpecies())

    for y in range(0, reaction.getNumProducts()):
      if reaction.getProduct(y).getSpecies() not in species_metab:
        species.add(reaction.getProduct(y).getSpecies())

    for y in range(0, reaction.getNumModifiers()):
      if reaction.getModifier(y).getSpecies() not in species_metab:
        species.add(reaction.getModifier(y).getSpecies())

    reactions.add(x) 

  #print (len(species_metab)) 
  #print (len(species_in_metab)) 

  # if metabolism is not included in the selected processes find only the common species
  # with the metabolism in order to initialize with higher value
  if include_Metabolism == False :
    species_metab = species.intersection(species_in_metab)

  #print (len(species_metab)) 
  assignment_rules = set()
  avogadro_con = 6.02214e23 
  volume_space = 5.23e-16
  print('<?xml version="1.0" encoding="UTF-8"?>')
  print('<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" level="3" version="1">')
  print('<model timeUnits="second" substanceUnits="item" extentUnits="item" lengthUnits="metre" areaUnits="metre2" volumeUnits="litre">')
  print('<listOfUnitDefinitions>')
  print('<unitDefinition id="rate">')
  print('<listOfUnits>')
  print('<unit kind="item" exponent="0" scale="0" multiplier="1" />')
  print('<unit kind="second" exponent="-1" scale="0" multiplier="1" />')
  print('</listOfUnits>')
  print('</unitDefinition>')
  print('<unitDefinition id="metre2">')
  print('<listOfUnits>')
  print('<unit kind="metre" exponent="2" scale="0" multiplier="1" />')
  print('</listOfUnits>')
  print('</unitDefinition>')
  print('</listOfUnitDefinitions>')
  print (model.getListOfCompartments().toSBML())
  print('<listOfSpecies>')
  
  #no input filename
  if input_used == False :
    #initialize and print metabolism species
    for x in species_metab:
      initvalue = numpy.random.randint(min_range_value, max_range_value + 1)
      model.getSpecies(x).setInitialAmount(float(initvalue))
      print ("  " + model.getSpecies(x).toSBML())
      if include_Metabolism == False :
        species.remove (x)

  #initialize and print species
  for x in species:
    # if the input for the metabolism species is defined
    if input_used == True :
      #print(str(x))
      if model_input.getSpecies(x) != None and x in species_metab:
        initvalue = model_input.getSpecies(x).getInitialAmount()*10
      else:
        if max_range_value > 500000:
          initvalue = numpy.random.randint(min_range_value/60, max_range_value/250 + 1)
        else:
          initvalue = numpy.random.randint(min_range_value/6, max_range_value/25 + 1)
    else:
      if max_range_value > 500000:
        initvalue = numpy.random.randint(min_range_value/60, max_range_value/250 + 1)
      else:
        initvalue = numpy.random.randint(min_range_value/6, max_range_value/25 + 1)
    model.getSpecies(x).setInitialAmount(float(initvalue))
    print ("  " + model.getSpecies(x).toSBML())

  print('</listOfSpecies>')
  print('<listOfReactions>')
  # add reaction rate formula
  keys = reactions_dict.keys()
  for key in keys:
    reaction_list = reactions_dict.get(key)
    
    #unique k for each cellular process
    mu = (min_k_value + max_k_value)/2
    sigma = (max_k_value - mu)/3 # +_3sigma (three standard deviations account for 99.73% probability)
    k = min_k_value -1
    while k < min_k_value or k > max_k_value:
      k = numpy.random.normal(mu, sigma)
    for rx in reaction_list:
      reaction = model.getReaction(rx)
      kl = reaction.createKineticLaw()


      # Delete H2O and H
      reactants = reaction.getListOfReactants()
      products = reaction.getListOfProducts()
      products.remove("c_H")
      products.remove("c_H2O")
      products.remove("e_H")
      products.remove("e_H2O")
      reactants.remove("c_H")
      reactants.remove("c_H2O")
      reactants.remove("e_H")
      reactants.remove("e_H2O")

      #print("key" + str(key))

      #protein degradation
      if key == 6:
        total = 0
        for y in range(0, reaction.getNumProducts()):
          species_name = reaction.getProduct(y).getSpecies()
          if (species_name[-4:] == "_ARG" or species_name[-4:] == "_HIS" 
          or species_name[-4:] == "_LYS" or species_name[-4:] == "_ASP" 
          or species_name[-4:] == "_GLU" or species_name[-4:] == "_SER" 
          or species_name[-4:] == "_THR" or species_name[-4:] == "_ASN" 
          or species_name[-4:] == "_GLN" or species_name[-4:] == "_CYS" 
          or species_name[-4:] == "_GLY" or species_name[-4:] == "_PRO" 
          or species_name[-4:] == "_ALA" or species_name[-4:] == "_VAL"
          or species_name[-4:] == "_ILE" or species_name[-4:] == "_LEU" 
          or species_name[-4:] == "_MET" or species_name[-4:] == "_PHE" 
          or species_name[-4:] == "_TYR" or species_name[-4:] == "_TRP" ):
            total=total+ reaction.getProduct(y).getStoichiometry()
        frac_deg=0.5/total
        
        k = frac_deg

        # Add k as a local parameter to reaction
        para = kl.createParameter()
        para.setId("k")
        para.setValue(k)
        para.setUnits("rate")
        para.setConstant(True)
        tempstring =""
        cnt_tps = 0
        min_conc = 100000000000.00
        name_conc = ""
        mathXMLString = """<math xmlns="http://www.w3.org/1998/Math/MathML">
        <apply>
          <times/>
            <ci> k </ci>"""
        # Add reactants
        for y in range(0, reaction.getNumReactants()):
          name_species=reaction.getReactant(y).getSpecies()
          species_conc = model.getSpecies(name_species).getInitialAmount()*1000 / (avogadro_con * volume_space)

          if (name_species != "c_GTP" and name_species != "c_ATP"):

            conc_string = """
              <apply>
               <divide/>
               <apply>
                 <times/>"""
            conc_string =conc_string + "<ci>" + reaction.getReactant(y).getSpecies() + "</ci>"
            conc_string =conc_string + """<cn type="integer"> 1000 </cn>
                </apply>
              <apply>
               <times/>
                <cn type="e-notation"> 6.02214 <sep/> 23 </cn>
                <cn type="e-notation"> 5.23 <sep/> -16 </cn>
              </apply>
             </apply>
            """

            if reaction.getReactant(y).getStoichiometry() > 10:
              print("Stoichiometry in reactants greater than 10 in reaction " + reaction.getIdAttribute() + " in species " + reaction.getReactant(y).getSpecies() )
              return 1
            elif reaction.getReactant(y).getStoichiometry() == 0:
              print("Stoichiometry 0 in reaction " + reaction.getIdAttribute() + " in species " + reaction.getReactant(y).getSpecies() )
              return 1
            if reaction.getReactant(y).getStoichiometry() == 1:
              mathXMLString =mathXMLString + conc_string
            else:
              mathXMLString =mathXMLString + "<apply><power/>" + conc_string + "<cn type=\"integer\">"+ str(reaction.getReactant(y).getStoichiometry())+"</cn></apply>"

          else:
            if species_conc < min_conc :
              min_conc = species_conc
              name_conc = reaction.getReactant(y).getSpecies()
            cnt_tps = cnt_tps + 1

        # Add modifiers
        for y in range(0, reaction.getNumModifiers()):
          name_species=reaction.getModifier(y).getSpecies()
          species_conc = model.getSpecies(name_species).getInitialAmount()*1000 / (avogadro_con * volume_space)
          conc_string = """
              <apply>
               <divide/>
               <apply>
                 <times/>"""
          conc_string =conc_string + "<ci>" + reaction.getModifier(y).getSpecies() + "</ci>"
          conc_string =conc_string + """<cn type="integer"> 1000 </cn>
                </apply>
              <apply>
               <times/>
                <cn type="e-notation"> 6.02214 <sep/> 23 </cn>
                <cn type="e-notation"> 5.23 <sep/> -16 </cn>
              </apply>
             </apply>
            """
          mathXMLString =mathXMLString  + conc_string


        if cnt_tps > 0 :

          conc_string = """
              <apply>
               <divide/>
               <apply>
                 <times/>"""
          conc_string =conc_string + "<ci>" + name_conc + "</ci>"
          conc_string =conc_string + """<cn type="integer"> 1000 </cn>
                </apply>
              <apply>
               <times/>
                <cn type="e-notation"> 6.02214 <sep/> 23 </cn>
                <cn type="e-notation"> 5.23 <sep/> -16 </cn>
              </apply>
             </apply>
            """

          mathXMLString =mathXMLString  + conc_string

        mathXMLString =mathXMLString + """
        </apply>
        </math>
        """

      
      # RNA degradation
      elif key == 4:
        total = 0
        #calculate k_RNA_deg
        for y in range(0, reaction.getNumProducts()):
          species_name = reaction.getProduct(y).getSpecies()
          if (species_name == "c_AMP" or species_name == "c_GMP" 
          or species_name == "c_CMP" or species_name == "c_UMP"):
            total=total+ reaction.getProduct(y).getStoichiometry()
        k_RNA_deg = 5/total
        #print ("\n total: " + str(total))
        # Add k_RNA_deg as a local parameter to reaction
        para = kl.createParameter()
        para.setId("k")
        para.setValue(k_RNA_deg)
        para.setUnits("rate")
        mathXMLString = """<math xmlns="http://www.w3.org/1998/Math/MathML">
        <apply>
          <times/>
            <ci> k </ci>"""

        for y in range(0, reaction.getNumReactants()):
          name_species=reaction.getReactant(y).getSpecies()
          species_conc = model.getSpecies(name_species).getInitialAmount()*1000 / (avogadro_con * volume_space)
          conc_string = """
              <apply>
               <divide/>
               <apply>
                 <times/>"""
          conc_string =conc_string + "<ci>" + reaction.getReactant(y).getSpecies() + "</ci>"
          conc_string =conc_string + """<cn type="integer"> 1000 </cn>
                </apply>
              <apply>
               <times/>
                <cn type="e-notation"> 6.02214 <sep/> 23 </cn>
                <cn type="e-notation"> 5.23 <sep/> -16 </cn>
              </apply>
             </apply>
            """

          if reaction.getReactant(y).getStoichiometry() > 1:
            print("Stoichiometry in reactants greater than 1 in reaction " + reaction.getIdAttribute() + " in species " + reaction.getReactant(y).getSpecies() )
            return 1
          elif reaction.getReactant(y).getStoichiometry() == 0:
            print("Stoichiometry 0 in reaction " + reaction.getIdAttribute() + " in species " + reaction.getReactant(y).getSpecies() )
            return 1
          if reaction.getReactant(y).getStoichiometry() == 1:
            mathXMLString =mathXMLString + conc_string
          else:
            mathXMLString =mathXMLString + "<apply><power/>" + conc_string + "<cn type=\"integer\">"+ str(reaction.getReactant(y).getStoichiometry())+"</cn></apply>"

        # Add modifiers
        for y in range(0, reaction.getNumModifiers()):
          name_species=reaction.getModifier(y).getSpecies()
          species_conc = model.getSpecies(name_species).getInitialAmount()*1000 / (avogadro_con * volume_space)
          conc_string = """
              <apply>
               <divide/>
               <apply>
                 <times/>"""
          conc_string =conc_string + "<ci>" + reaction.getModifier(y).getSpecies() + "</ci>"
          conc_string =conc_string + """<cn type="integer"> 1000 </cn>
                </apply>
              <apply>
               <times/>
                <cn type="e-notation"> 6.02214 <sep/> 23 </cn>
                <cn type="e-notation"> 5.23 <sep/> -16 </cn>
              </apply>
             </apply>
            """
          mathXMLString =mathXMLString  + conc_string

        mathXMLString =mathXMLString + """
        </apply>
        </math>
        """






      # the rest of processes
      else:
        total = 0
        #calculate k_RNA_deg
        for y in range(0, reaction.getNumReactants()):
          total=total+ reaction.getReactant(y).getStoichiometry()
        #print(total)
        # first approach sum of stoichiometries < 7
        if (total < max_sum_stoich ):
          if (approach1 == True): 
            # Add k as a local parameter to reaction
            para = kl.createParameter()
            para.setId("k")
            para.setValue(k)
            para.setUnits("rate")
            para.setConstant(True)
            tempstring =""
            cnt_tps = 0
            min_conc = 100000000000.00
            name_conc = ""
            total_conc_stoic = 0
            mathXMLString = """<math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
                <ci> k </ci>"""
            for y in range(0, reaction.getNumReactants()):
              name_species=reaction.getReactant(y).getSpecies()
              species_conc = model.getSpecies(name_species).getInitialAmount()*1000 / (avogadro_con * volume_space)

              conc_string = """
                  <apply>
                  <divide/>
                  <apply>
                    <times/>"""
              conc_string =conc_string + "<ci>" + reaction.getReactant(y).getSpecies() + "</ci>"
              conc_string =conc_string + """<cn type="integer"> 1000 </cn>
                    </apply>
                  <apply>
                  <times/>
                    <cn type="e-notation"> 6.02214 <sep/> 23 </cn>
                    <cn type="e-notation"> 5.23 <sep/> -16 </cn>
                  </apply>
                </apply>
                """
              if reaction.getReactant(y).getStoichiometry() == 1:
                mathXMLString =mathXMLString + conc_string
              else:
                mathXMLString =mathXMLString + "<apply><power/>" + conc_string + "<cn type=\"integer\">"+ str(reaction.getReactant(y).getStoichiometry())+"</cn></apply>"

            # Add modifiers
            for y in range(0, reaction.getNumModifiers()):
              name_species=reaction.getModifier(y).getSpecies()
              species_conc = model.getSpecies(name_species).getInitialAmount()*1000 / (avogadro_con * volume_space)
              conc_string = """
                  <apply>
                  <divide/>
                  <apply>
                    <times/>"""
              conc_string =conc_string + "<ci>" + reaction.getModifier(y).getSpecies() + "</ci>"
              conc_string =conc_string + """<cn type="integer"> 1000 </cn>
                    </apply>
                  <apply>
                  <times/>
                    <cn type="e-notation"> 6.02214 <sep/> 23 </cn>
                    <cn type="e-notation"> 5.23 <sep/> -16 </cn>
                  </apply>
                </apply>
                """
              mathXMLString =mathXMLString  + conc_string

            mathXMLString =mathXMLString + """
            </apply>
            </math>
            """
        # three approach sum of stoichiometries >= 7  
        elif (total >= max_sum_stoich ):
          if (approach3 == True):
            #Replication
            if key == 1:
              lp=0
              k_pol=0
              for y in range(0, reaction.getNumReactants()):
                species_name = reaction.getReactant(y).getSpecies()
                if (species_name == "c_DCTP"  or species_name == "c_DATP"
                or species_name == "c_DGTP" or species_name == "c_DTTP"
                or species_name== "c_DUTP" or species_name[-8:] == "_MONOMER" 
                or species_name[-8:] == "_HEXAMER" or species_name[-9:] == "_TETRAMER"):
                  lp=lp+ reaction.getReactant(y).getStoichiometry()
              if lp == 0:
                print("lp 0 in reaction " + reaction.getIdAttribute() )
                return 1    
              k_pol = 1000/lp
            #Transcription  
            elif key == 2:
              lp=0
              k_pol=0
              for y in range(0, reaction.getNumReactants()):
                species_name = reaction.getReactant(y).getSpecies()
                if (species_name == "c_GTP" or species_name == "c_ATP" 
                or species_name == "c_CTP"  or species_name == "c_UTP"
                or species_name == "c_TTP"):
                  lp=lp+ reaction.getReactant(y).getStoichiometry()
              k_pol = 60/lp
            #Translation  
            elif key == 5:
              lp=0
              k_pol=0 
              for y in range(0, reaction.getNumReactants()):
                species_name = reaction.getReactant(y).getSpecies()
                if (species_name[-4:] == "_ARG" or species_name[-4:] == "_HIS" 
                or species_name[-4:] == "_LYS" or species_name[-4:] == "_ASP" 
                or species_name[-4:] == "_GLU" or species_name[-4:] == "_SER" 
                or species_name[-4:] == "_THR" or species_name[-4:] == "_ASN" 
                or species_name[-4:] == "_GLN" or species_name[-4:] == "_CYS" 
                or species_name[-4:] == "_GLY" or species_name[-4:] == "_PRO" 
                or species_name[-4:] == "_ALA" or species_name[-4:] == "_VAL"
                or species_name[-4:] == "_ILE" or species_name[-4:] == "_LEU" 
                or species_name[-4:] == "_MET" or species_name[-4:] == "_PHE" 
                or species_name[-4:] == "_TYR" or species_name[-4:] == "_TRP" ):
                  lp = lp + reaction.getReactant(y).getStoichiometry()
              k_pol = 20/lp
            #Protein maturation  
            elif key == 8:
              lp=0
              k_pol=0
              for y in range(0, reaction.getNumReactants()):
                species_name = reaction.getReactant(y).getSpecies()
                if (species_name == "c_GTP" or species_name == "c_ATP" 
                or species_name == "c_CTP"  or species_name == "c_UTP"
                or species_name == "c_TTP"):
                  lp=lp+ reaction.getReactant(y).getStoichiometry()
              k_pol = 0.05/lp
            #Macromolecular complexation  
            elif key == 9:
              lp=0
              k_pol=0
              for y in range(0, reaction.getNumReactants()):
                species_name = reaction.getReactant(y).getSpecies()
                if (species_name[-8:] == ("_MONOMER")):
                  lp=lp+ reaction.getReactant(y).getStoichiometry()
              k_pol = 20/lp
            #RNA processing  
            elif key == 13:
              lp=0
              k_pol=0
              for y in range(0, reaction.getNumReactants()):
                species_name = reaction.getReactant(y).getSpecies()
                if (species_name == "c_GTP" or species_name == "c_ATP" 
                or species_name == "c_CTP"  or species_name == "c_UTP"
                or species_name == "c_TTP"):
                  lp=lp+ reaction.getReactant(y).getStoichiometry()
              k_pol = 0.05/lp
            else:
              print("More cellular processes with reactions with sum of stoichiometries >= " + max_sum_stoich + " in "  + reaction.getIdAttribute()  )
              return 1

            k = k_pol
            # Add k as a local parameter to reaction
            para = kl.createParameter()
            para.setId("k")
            para.setValue(k)
            para.setUnits("rate")
            para.setConstant(True)
            tempstring =""
            cnt_tps = 0
            cnt_tpsA = 0
            min_conc = 100000000000.00
            min_concA = 100000000000.00
            name_conc = ""
            name_concA = ""
            mathXMLString = """<math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply>
              <times/>
                <ci> k </ci>"""
            for y in range(0, reaction.getNumReactants()):
              name_species=reaction.getReactant(y).getSpecies()
              species_conc = model.getSpecies(name_species).getInitialAmount()*1000 / (avogadro_con * volume_space)
              #print (name_species)
              
              #Replication
              if key == 1:
                if (name_species != "c_DCTP" and name_species != "c_DATP"
                and name_species != "c_DGTP" and name_species != "c_DTTP"
                and name_species != "c_DUTP" and name_species != "c_ATP"):
              
                  conc_string = """
                    <apply>
                    <divide/>
                    <apply>
                      <times/>"""
                  conc_string =conc_string + "<ci>" + reaction.getReactant(y).getSpecies() + "</ci>"
                  conc_string =conc_string + """<cn type="integer"> 1000 </cn>
                      </apply>
                    <apply>
                    <times/>
                      <cn type="e-notation"> 6.02214 <sep/> 23 </cn>
                      <cn type="e-notation"> 5.23 <sep/> -16 </cn>
                    </apply>
                  </apply>
                  """

                  if reaction.getReactant(y).getStoichiometry() > 8:
                    print("Stoichiometry in reactants greater than 8 in reaction " + reaction.getIdAttribute() + " in species " + reaction.getReactant(y).getSpecies() )
                    return 1
                  if reaction.getReactant(y).getStoichiometry() == 1:
                    mathXMLString =mathXMLString + conc_string
                  else:
                    mathXMLString =mathXMLString + "<apply><power/>" + conc_string + "<cn type=\"integer\">"+ str(reaction.getReactant(y).getStoichiometry())+"</cn></apply>"

                else:
                  if species_conc < min_conc :
                    min_conc = species_conc
                    name_conc = reaction.getReactant(y).getSpecies()
                  cnt_tps = cnt_tps + 1
              #Translation  
              elif key == 5: 

                if (name_species[-4:] != "_ARG" and name_species[-4:] != "_HIS" 
                and name_species[-4:] != "_LYS" and name_species[-4:] != "_ASP" 
                and name_species[-4:] != "_GLU" and name_species[-4:] != "_SER" 
                and name_species[-4:] != "_THR" and name_species[-4:] != "_ASN" 
                and name_species[-4:] != "_GLN" and name_species[-4:] != "_CYS" 
                and name_species[-4:] != "_GLY" and name_species[-4:] != "_PRO" 
                and name_species[-4:] != "_ALA" and name_species[-4:] != "_VAL"
                and name_species[-4:] != "_ILE" and name_species[-4:] != "_LEU" 
                and name_species[-4:] != "_MET" and name_species[-4:] != "_PHE" 
                and name_species[-4:] != "_TYR" and name_species[-4:] != "_TRP"
                and name_species != "c_CTP" and name_species != "c_ATP"
                and name_species != "c_GTP" and name_species != "c_TTP"
                and name_species!= "c_UTP" ):
                  
                  conc_string = """
                    <apply>
                    <divide/>
                    <apply>
                      <times/>"""
                  conc_string =conc_string + "<ci>" + reaction.getReactant(y).getSpecies() + "</ci>"
                  conc_string =conc_string + """<cn type="integer"> 1000 </cn>
                      </apply>
                    <apply>
                    <times/>
                      <cn type="e-notation"> 6.02214 <sep/> 23 </cn>
                      <cn type="e-notation"> 5.23 <sep/> -16 </cn>
                    </apply>
                  </apply>
                  """

                  if reaction.getReactant(y).getStoichiometry() > 4:
                    print("Stoichiometry in reactants greater than 4 in reaction " + reaction.getIdAttribute() + " in species " + reaction.getReactant(y).getSpecies() )
                    return 1
                  if reaction.getReactant(y).getStoichiometry() == 1:
                    mathXMLString =mathXMLString + conc_string
                  else:
                    mathXMLString =mathXMLString + "<apply><power/>" + conc_string + "<cn type=\"integer\">"+ str(reaction.getReactant(y).getStoichiometry())+"</cn></apply>"

                else:
                  if (name_species[-4:] == "_ARG" or name_species[-4:] == "_HIS" 
                  or name_species[-4:] == "_LYS" or name_species[-4:] == "_ASP" 
                  or name_species[-4:] == "_GLU" or name_species[-4:] == "_SER" 
                  or name_species[-4:] == "_THR" or name_species[-4:] == "_ASN" 
                  or name_species[-4:] == "_GLN" or name_species[-4:] == "_CYS" 
                  or name_species[-4:] == "_GLY" or name_species[-4:] == "_PRO" 
                  or name_species[-4:] == "_ALA" or name_species[-4:] == "_VAL"
                  or name_species[-4:] == "_ILE" or name_species[-4:] == "_LEU" 
                  or name_species[-4:] == "_MET" or name_species[-4:] == "_PHE" 
                  or name_species[-4:] == "_TYR" or name_species[-4:] == "_TRP" ):
                    if species_conc < min_concA :
                      min_concA = species_conc
                      name_concA = reaction.getReactant(y).getSpecies()
                    cnt_tpsA = cnt_tpsA + 1 

                  elif (name_species == "c_GTP" or name_species == "c_ATP" 
                    or name_species == "c_CTP"  or name_species == "c_UTP"
                    or name_species == "c_TTP"): 
                    
                    if species_conc < min_conc :
                      min_conc = species_conc
                      name_conc = reaction.getReactant(y).getSpecies()
                    cnt_tps = cnt_tps + 1
                    
                  

              #Macromolecular complexation  
              elif key == 9:
                if (name_species[-8:] != "_MONOMER" and name_species != "c_CTP" 
                and name_species != "c_ATP" and name_species != "c_GTP" 
                and name_species != "c_TTP" and name_species!= "c_UTP"): 
                  conc_string = """
                    <apply>
                    <divide/>
                    <apply>
                      <times/>"""
                  conc_string =conc_string + "<ci>" + reaction.getReactant(y).getSpecies() + "</ci>"
                  conc_string =conc_string + """<cn type="integer"> 1000 </cn>
                      </apply>
                    <apply>
                    <times/>
                      <cn type="e-notation"> 6.02214 <sep/> 23 </cn>
                      <cn type="e-notation"> 5.23 <sep/> -16 </cn>
                    </apply>
                  </apply>
                  """
                  if reaction.getReactant(y).getStoichiometry() > 4:
                    print("Stoichiometry in reactants greater than 4 in reaction " + reaction.getIdAttribute() + " in species " + reaction.getReactant(y).getSpecies() )
                    return 1
                  if reaction.getReactant(y).getStoichiometry() == 1:
                    mathXMLString =mathXMLString + conc_string
                  else:
                    mathXMLString =mathXMLString + "<apply><power/>" + conc_string + "<cn type=\"integer\">"+ str(reaction.getReactant(y).getStoichiometry())+"</cn></apply>"

                else:
                  if (name_species[-8:] == "_MONOMER" ): 
                    if species_conc < min_conc :
                      min_conc = species_conc
                      name_conc = reaction.getReactant(y).getSpecies()
                    cnt_tps = cnt_tps + 1 

                  elif (name_species == "c_GTP" or name_species == "c_ATP" 
                    or name_species == "c_CTP"  or name_species == "c_UTP"
                    or name_species == "c_TTP"):
                    if species_conc < min_concA :
                      min_concA = species_conc
                      name_concA = reaction.getReactant(y).getSpecies()
                    cnt_tpsA = cnt_tpsA + 1  
                  

              #Transcription, Protein maturation  and RNA processing 
              else:
                if (name_species != "c_CTP" and name_species != "c_ATP"
                and name_species != "c_GTP" and name_species != "c_TTP"
                and name_species!= "c_UTP"):
              
                  conc_string = """
                    <apply>
                    <divide/>
                    <apply>
                      <times/>"""
                  conc_string =conc_string + "<ci>" + reaction.getReactant(y).getSpecies() + "</ci>"
                  conc_string =conc_string + """<cn type="integer"> 1000 </cn>
                      </apply>
                    <apply>
                    <times/>
                      <cn type="e-notation"> 6.02214 <sep/> 23 </cn>
                      <cn type="e-notation"> 5.23 <sep/> -16 </cn>
                    </apply>
                  </apply>
                  """

                  if reaction.getReactant(y).getStoichiometry() > 4:
                    print("Stoichiometry in reactants greater than 4 in reaction " + reaction.getIdAttribute() + " in species " + reaction.getReactant(y).getSpecies() )
                    return 1
                  if reaction.getReactant(y).getStoichiometry() == 1:
                    mathXMLString =mathXMLString + conc_string
                  else:
                    mathXMLString =mathXMLString + "<apply><power/>" + conc_string + "<cn type=\"integer\">"+ str(reaction.getReactant(y).getStoichiometry())+"</cn></apply>"

                else:
                  if species_conc < min_conc :
                    min_conc = species_conc
                    name_conc = reaction.getReactant(y).getSpecies()
                  cnt_tps = cnt_tps + 1

            # Add modifiers
            for y in range(0, reaction.getNumModifiers()):
              name_species=reaction.getModifier(y).getSpecies()
              species_conc = model.getSpecies(name_species).getInitialAmount()*1000 / (avogadro_con * volume_space)
              conc_string = """
                  <apply>
                  <divide/>
                  <apply>
                    <times/>"""
              conc_string =conc_string + "<ci>" + reaction.getModifier(y).getSpecies() + "</ci>"
              conc_string =conc_string + """<cn type="integer"> 1000 </cn>
                    </apply>
                  <apply>
                  <times/>
                    <cn type="e-notation"> 6.02214 <sep/> 23 </cn>
                    <cn type="e-notation"> 5.23 <sep/> -16 </cn>
                  </apply>
                </apply>
                """
              mathXMLString =mathXMLString  + conc_string

            if cnt_tps > 0 :
              conc_string = """
                  <apply>
                  <divide/>
                  <apply>
                    <times/>"""
              conc_string =conc_string + "<ci>" + name_conc + "</ci>"
              conc_string =conc_string + """<cn type="integer"> 1000 </cn>
                    </apply>
                  <apply>
                  <times/>
                    <cn type="e-notation"> 6.02214 <sep/> 23 </cn>
                    <cn type="e-notation"> 5.23 <sep/> -16 </cn>
                  </apply>
                </apply>
                """
              mathXMLString =mathXMLString + conc_string

            if cnt_tpsA > 0 :
              conc_string = """
                  <apply>
                  <divide/>
                  <apply>
                    <times/>"""
              conc_string =conc_string + "<ci>" + name_concA + "</ci>"
              conc_string =conc_string + """<cn type="integer"> 1000 </cn>
                    </apply>
                  <apply>
                  <times/>
                    <cn type="e-notation"> 6.02214 <sep/> 23 </cn>
                    <cn type="e-notation"> 5.23 <sep/> -16 </cn>
                  </apply>
                </apply>
                """
              mathXMLString =mathXMLString + conc_string  



            mathXMLString =mathXMLString + """
            </apply>
            </math>
            """


      astMath = readMathMLFromString(mathXMLString)
      kl.setMath(astMath)
      if reaction.getNumReactants() != 0 and reaction.getNumProducts()!=0:
        print (reaction.toSBML())
  print('</listOfReactions>')
  if model.getListOfRules().size() > 0 :
    print (model.getListOfRules().toSBML())
  if model.getListOfParameters().size() > 0 :
    print (model.getListOfParameters().toSBML())
  print('</model>')
  print('</sbml>')

  if args.output_filename != None :
    #f.write("Woops! I have deleted the content!")
    sys.stdout = original_stdout # Reset the standard output to its original value
    f.close()

  return 0


if __name__ == '__main__':
  main(sys.argv)
