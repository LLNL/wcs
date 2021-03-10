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
  optionalNamed.add_argument('-o', metavar='<output_filename>', nargs=1, help='Specify the output filename', dest="output_filename")
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
species  (default from 10 to 1000)''',
  dest="range_values", default=["10", "1000"])
  optionalNamed.add_argument('-n', metavar=('<min>', '<max>'), nargs=2, type=float, help='''Specify the range of the normal distribution 
for the initialization of reaction rate coefficient k used in the
reaction rate formulas (default from 2.5 to 7.5)''',
  dest="normal_distr", default=["2.5", "7.5"])

  original_stdout = sys.stdout # Save a reference to the original standard output

  try:
    args = parser.parse_args()

    #print('Range:', args.range_values)
    if args.range_values != None :
      if (args.range_values[0] > args.range_values[1] or args.range_values[0] < 0 ) :
        print('Error: the first argument of -r has to be greater than 0 and smaller than the second argument')
        sys.exit()
      min_range_value = args.range_values[0]
      max_range_value = args.range_values[1]
      #print ("\n" + str(min_range_value) + " " + str(max_range_value))

    if args.normal_distr != None :
      if (args.normal_distr[0] > args.normal_distr[1] or args.normal_distr[0] < 0 ) :
        print('Error: the first argument of -n has to be greater than 0 and smaller than the second argument')
        sys.exit()
      min_k_value = args.normal_distr[0]
      max_k_value = args.normal_distr[1]
    #print ("\n" + str(min_k_value) + " " + str(max_k_value))

    processes = set ()
    #print('Processes:', args.proc)
    if args.proc != None :
      for opt in args.proc:
        if (int(opt) > 15 or int(opt) < 0 ) :
          print('Error: the numbers corresponding to processes have to be between 1 to 15. Run \"python generate_network.py -h\" to see the number for each process.')
          sys.exit()
        processes.add(opt)

    #for opt in processes:
    #    print ("\n" + str(opt))

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

  reactions_ids = set ()
  reactions_dict = {}
  if len(processes) > 0 :
    processes_names = ["Replication", "Transcription","Protein-DNA Interaction ","RNA Degradation","Translation"
    ,"Protein Degradation","Metabolism","Protein Maturation","Macromolecular Complexation","Protein-DNA Interaction"
    ,"tRNA Aminoacylation","Terminal Organelle Assembly","RNA Processing","Replication Initiation","Cellular Division"]

  #print('\n Type ='+ str(len(processes)))
  for elem in root:
    for subelem in elem:
      for subsubelem in subelem:
        for  proc in processes:
          if (subsubelem.text == processes_names[int(proc)-1]):
            reactions_ids.add(str(subelem.attrib.get('id')))
            reactions_dict.setdefault(int(proc), [])
            reactions_dict[int(proc)].append(str(subelem.attrib.get('id')))

  #print(" Reactions:" + str(len(reactions_dict)))
  #print(" Reactions:" + str(len(reactions_ids)))
  #print(reactions_dict)

  #keys = reactions_dict.keys()
  #for x in keys:
    #print ("\n key " + str(x) + " :")
    #rlist = reactions_dict.get(x)
    #for y in rlist:
    #  print (str(y) + ", ")

  if document.getNumErrors() > 0:
    print("Encountered the following SBML errors:" )
    document.printErrors()
    return 1

  model = document.getModel()

  if model is None:
    print("No model present." )
    return 1

  species = set()
  reactions = set ()
  for x in reactions_ids:
    reaction = model.getReaction(x)
    for y in range(0, reaction.getNumReactants()):
      species.add(reaction.getReactant(y).getSpecies())
    for y in range(0, reaction.getNumProducts()):
      species.add(reaction.getProduct(y).getSpecies())
    for y in range(0, reaction.getNumModifiers()):
      species.add(reaction.getModifier(y).getSpecies())
    reactions.add(x)
  assignment_rules = set()
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

  # find maximum reactant stoichiometry
  # max_stoich = 0
  # keys = reactions_dict.keys()
  # for key in keys:
  #   reaction_list = reactions_dict.get(key)
  #   for rx in reaction_list:
  #     reaction = model.getReaction(rx)
  #     for y in range(0, reaction.getNumReactants()):
  #       if (reaction.getReactant(y).getSpecies() != "c_H2O" and reaction.getReactant(y).getSpecies() != "c_H"):
  #         if reaction.getReactant(y).getStoichiometry() > max_stoich:
  #           #print('Species: ' + reaction.getReactant(y).getSpecies()+  " " + str(reaction.getReactant(y).getStoichiometry()))
  #           max_stoich = reaction.getReactant(y).getStoichiometry()
  #   #print('Reactions: ' + str(len(reaction_list)))
  # #print('Max stoich: '+ str(max_stoich))

  #initialize and print species
  for x in species:
    initvalue = numpy.random.randint(min_range_value, max_range_value + 1)
    #initvalue = numpy.random.randint(min_range_value * max_stoich, max_range_value * max_stoich + 1)
    model.getSpecies(x).setInitialAmount(float(initvalue))
    print ("  " + model.getSpecies(x).toSBML())
  print('</listOfSpecies>')
  print('<listOfReactions>')
  # add reaction rate formula
  keys = reactions_dict.keys()
  for key in keys:
    reaction_list = reactions_dict.get(key)
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
          total=total+ reaction.getProduct(y).getStoichiometry()
        frac_deg=30/total
        #k=2
        mu = (min_k_value + max_k_value)/2
        sigma = (max_k_value - mu)/3 # +_3sigma (three standard deviations account for 99.73% probability)
        k_in = min_k_value -1
        while k_in < min_k_value or k_in > max_k_value:
          k_in = numpy.random.normal(mu, sigma)
          #print ("\n" + str(k_in))

        k = k_in * frac_deg

        # Add k as a local parameter to reaction
        para = kl.createParameter()
        para.setId("k")
        para.setValue(k)
        para.setUnits("rate")
        para.setConstant(True)
        tempstring =""
        cnt_tps = 0
        min_conc = 100000000000
        name_conc = ""
        total_conc_stoic = 0
        mathXMLString = """<math xmlns="http://www.w3.org/1998/Math/MathML">
        <apply>
          <times/>
            <ci> k </ci>
            <apply>
            <plus/>"""
        # Add reactants
        for y in range(0, reaction.getNumReactants()):
          name_species=reaction.getReactant(y).getSpecies()
          species_conc = model.getSpecies(name_species).getInitialAmount()*1000 / (6e23*5.23e-16)

          if (reaction.getReactant(y).getSpecies() != "c_GTP" and reaction.getReactant(y).getSpecies() != "c_ATP"
          and reaction.getReactant(y).getSpecies() != "c_CTP" and reaction.getReactant(y).getSpecies() != "c_UTP"
          and reaction.getReactant(y).getSpecies() != "c_DCTP" and reaction.getReactant(y).getSpecies() != "c_DATP"
          and reaction.getReactant(y).getSpecies() != "c_DGTP" and reaction.getReactant(y).getSpecies() != "c_DTTP"
          and reaction.getReactant(y).getSpecies() != "c_TTP" and reaction.getReactant(y).getSpecies() != "c_UTP"):

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
                <cn type="e-notation"> 6 <sep/> 23 </cn>
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
            total_conc_stoic =  total_conc_stoic + reaction.getReactant(y).getStoichiometry()
            cnt_tps = cnt_tps + 1

        # Add modifiers
        for y in range(0, reaction.getNumModifiers()):
          name_species=reaction.getModifier(y).getSpecies()
          species_conc = model.getSpecies(name_species).getInitialAmount()*1000 / (6e23*5.23e-16)
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
                <cn type="e-notation"> 6 <sep/> 23 </cn>
                <cn type="e-notation"> 5.23 <sep/> -16 </cn>
              </apply>
             </apply>
            """
          mathXMLString =mathXMLString  + conc_string

        # close the total of the reactants and modifiers
        mathXMLString =mathXMLString +"""
             </apply>
            """

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
                <cn type="e-notation"> 6 <sep/> 23 </cn>
                <cn type="e-notation"> 5.23 <sep/> -16 </cn>
              </apply>
             </apply>
            """

          if total_conc_stoic == 1:
            mathXMLString =mathXMLString  + conc_string
          else:
            mathXMLString =mathXMLString + "<apply><divide/>" +  conc_string + "<cn type=\"integer\">"+ str(total_conc_stoic)+"</cn></apply>"



        mathXMLString =mathXMLString + """
        </apply>
        </math>
        """

      #translation
      elif key == 5:
        mu = (min_k_value + max_k_value)/2
        sigma = (max_k_value - mu)/3 # +_3sigma (three standard deviations account for 99.73% probability)
        k = min_k_value -1
        while k < min_k_value or k > max_k_value:
          k = numpy.random.normal(mu, sigma)

        # Add k as a local parameter to reaction
        para = kl.createParameter()
        para.setId("k")
        para.setValue(k)
        para.setUnits("rate")
        para.setConstant(True)
        tempstring =""
        cnt_tps = 0
        min_conc = 100000000000
        name_conc = ""
        total_conc_stoic = 0
        mathXMLString = """<math xmlns="http://www.w3.org/1998/Math/MathML">
        <apply>
          <times/>
            <ci> k </ci>"""
        # simple reaction
        if reaction.getNumReactants() < 5 :
          mathXMLString = mathXMLString + """
            <apply>
            <plus/>"""
          # for reactants
          for y in range(0, reaction.getNumReactants()):
            name_species=reaction.getReactant(y).getSpecies()
            species_conc = model.getSpecies(name_species).getInitialAmount()*1000 / (6e23*5.23e-16)

            if (reaction.getReactant(y).getSpecies() != "c_GTP" and reaction.getReactant(y).getSpecies() != "c_ATP"
            and reaction.getReactant(y).getSpecies() != "c_CTP" and reaction.getReactant(y).getSpecies() != "c_UTP"
            and reaction.getReactant(y).getSpecies() != "c_DCTP" and reaction.getReactant(y).getSpecies() != "c_DATP"
            and reaction.getReactant(y).getSpecies() != "c_DGTP" and reaction.getReactant(y).getSpecies() != "c_DTTP"
            and reaction.getReactant(y).getSpecies() != "c_TTP" and reaction.getReactant(y).getSpecies() != "c_UTP"):

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
                <cn type="e-notation"> 6 <sep/> 23 </cn>
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
              total_conc_stoic =  total_conc_stoic + reaction.getReactant(y).getStoichiometry()
              cnt_tps = cnt_tps + 1

          #For modifiers
          for y in range(0, reaction.getNumModifiers()):
            name_species=reaction.getModifier(y).getSpecies()
            species_conc = model.getSpecies(name_species).getInitialAmount()*1000 / (6e23*5.23e-16)

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
                <cn type="e-notation"> 6 <sep/> 23 </cn>
                <cn type="e-notation"> 5.23 <sep/> -16 </cn>
              </apply>
             </apply>
            """
            mathXMLString =mathXMLString + conc_string

          # close the total of the reactants and modifiers
          mathXMLString =mathXMLString +"""
             </apply>
            """

        # complex reaction
        else:
          for y in range(0, reaction.getNumReactants()):
            name_species=reaction.getReactant(y).getSpecies()
            species_conc = model.getSpecies(name_species).getInitialAmount()*1000 / (6e23*5.23e-16)
            if species_conc < min_conc :
              min_conc = species_conc
              name_conc = reaction.getReactant(y).getSpecies()
            total_conc_stoic =  total_conc_stoic + reaction.getReactant(y).getStoichiometry()
            cnt_tps = cnt_tps + 1

          # For modifiers
          for y in range(0, reaction.getNumModifiers()):
            name_species=reaction.getModifier(y).getSpecies()
            species_conc = model.getSpecies(name_species).getInitialAmount()*1000 / (6e23*5.23e-16)
            if species_conc < min_conc :
              min_conc = species_conc
              name_conc = reaction.getModifier(y).getSpecies()
            total_conc_stoic =  total_conc_stoic + 1
            cnt_tps = cnt_tps + 1

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
                <cn type="e-notation"> 6 <sep/> 23 </cn>
                <cn type="e-notation"> 5.23 <sep/> -16 </cn>
              </apply>
             </apply>
            """


          if total_conc_stoic == 1:
            mathXMLString =mathXMLString + conc_string
          else:
            mathXMLString =mathXMLString + "<apply><divide/>" + conc_string + "<cn type=\"integer\">"+ str(total_conc_stoic)+"</cn></apply>"

        mathXMLString =mathXMLString + """
        </apply>
        </math>
        """


      # RNA degradation
      elif key == 4:
        mu = (min_k_value + max_k_value)/2
        sigma = (max_k_value - mu)/3 # +_3sigma (three standard deviations account for 99.73% probability)

        ka = min_k_value -1
        while ka < min_k_value or ka > max_k_value:
          ka = numpy.random.normal(mu, sigma)

        kg = min_k_value -1
        while kg < min_k_value or kg > max_k_value:
          kg = numpy.random.normal(mu, sigma)

        kc = min_k_value -1
        while kc < min_k_value or kc > max_k_value:
          kc = numpy.random.normal(mu, sigma)

        ku = min_k_value -1
        while ku < min_k_value or ku > max_k_value:
          ku = numpy.random.normal(mu, sigma)

        total = 0
        for y in range(0, reaction.getNumProducts()):
          total=total+ reaction.getProduct(y).getStoichiometry()
        frac_deg=70/total
        #calculate k_RNA_deg
        k_RNA_deg_temp=0.0
        for y in range(0, reaction.getNumProducts()):
          if (reaction.getProduct(y).getSpecies() == "c_AMP"):
            k_RNA_deg_temp=k_RNA_deg_temp + (ka * (reaction.getProduct(y).getStoichiometry()/total))
          elif (reaction.getProduct(y).getSpecies() == "c_GMP"):
            k_RNA_deg_temp=k_RNA_deg_temp + (kg * (reaction.getProduct(y).getStoichiometry()/total))
          elif (reaction.getProduct(y).getSpecies() == "c_CMP"):
            k_RNA_deg_temp=k_RNA_deg_temp + (kc * (reaction.getProduct(y).getStoichiometry()/total))
          elif (reaction.getProduct(y).getSpecies() == "c_UMP"):
            k_RNA_deg_temp=k_RNA_deg_temp + (ku * (reaction.getProduct(y).getStoichiometry()/total))
        k_RNA_deg = frac_deg * k_RNA_deg_temp
        # Add k_RNA_deg as a local parameter to reaction
        para = kl.createParameter()
        para.setId("k")
        para.setValue(k_RNA_deg)
        para.setUnits("rate")
        mathXMLString = """<math xmlns="http://www.w3.org/1998/Math/MathML">
        <apply>
          <times/>
            <ci> k </ci>
            <apply>
            <plus/>"""
        for y in range(0, reaction.getNumReactants()):
          name_species=reaction.getReactant(y).getSpecies()
          species_conc = model.getSpecies(name_species).getInitialAmount()*1000 / (6e23*5.23e-16)
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
                <cn type="e-notation"> 6 <sep/> 23 </cn>
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
          species_conc = model.getSpecies(name_species).getInitialAmount()*1000 / (6e23*5.23e-16)
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
                <cn type="e-notation"> 6 <sep/> 23 </cn>
                <cn type="e-notation"> 5.23 <sep/> -16 </cn>
              </apply>
             </apply>
            """
          mathXMLString =mathXMLString  + conc_string

        # close the total of the reactants and modifiers
        mathXMLString =mathXMLString +"""
             </apply>
            """

        mathXMLString =mathXMLString + """
        </apply>
        </math>
        """






      # the rest of processes
      else:
        mu = (min_k_value + max_k_value)/2
        sigma = (max_k_value - mu)/3 # +_3sigma (three standard deviations account for 99.73% probability)
        k = min_k_value -1
        while k < min_k_value or k > max_k_value:
          k = numpy.random.normal(mu, sigma)

        # Add k as a local parameter to reaction
        para = kl.createParameter()
        para.setId("k")
        para.setValue(k)
        para.setUnits("rate")
        para.setConstant(True)
        tempstring =""
        cnt_tps = 0
        min_conc = 100000000000
        name_conc = ""
        total_conc_stoic = 0
        mathXMLString = """<math xmlns="http://www.w3.org/1998/Math/MathML">
         <apply>
           <times/>
             <ci> k </ci>
            <apply>
            <plus/>"""
        for y in range(0, reaction.getNumReactants()):
          name_species=reaction.getReactant(y).getSpecies()
          species_conc = model.getSpecies(name_species).getInitialAmount()*1000 / (6e23*5.23e-16)

          if (reaction.getReactant(y).getSpecies() != "c_GTP" and reaction.getReactant(y).getSpecies() != "c_ATP"
           and reaction.getReactant(y).getSpecies() != "c_CTP" and reaction.getReactant(y).getSpecies() != "c_UTP"
           and reaction.getReactant(y).getSpecies() != "c_DCTP" and reaction.getReactant(y).getSpecies() != "c_DATP"
           and reaction.getReactant(y).getSpecies() != "c_DGTP" and reaction.getReactant(y).getSpecies() != "c_DTTP"
           and reaction.getReactant(y).getSpecies() != "c_TTP" and reaction.getReactant(y).getSpecies() != "c_UTP"):
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
                <cn type="e-notation"> 6 <sep/> 23 </cn>
                <cn type="e-notation"> 5.23 <sep/> -16 </cn>
              </apply>
             </apply>
            """

            if reaction.getReactant(y).getStoichiometry() > 15:
              print("Stoichiometry in reactants greater than 15 in reaction " + reaction.getIdAttribute() + " in species " + reaction.getReactant(y).getSpecies() )
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
            total_conc_stoic =  total_conc_stoic + reaction.getReactant(y).getStoichiometry()
            cnt_tps = cnt_tps + 1

        # Add modifiers
        for y in range(0, reaction.getNumModifiers()):
          name_species=reaction.getModifier(y).getSpecies()
          species_conc = model.getSpecies(name_species).getInitialAmount()*1000 / (6e23*5.23e-16)
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
                <cn type="e-notation"> 6 <sep/> 23 </cn>
                <cn type="e-notation"> 5.23 <sep/> -16 </cn>
              </apply>
             </apply>
            """
          mathXMLString =mathXMLString  + conc_string


        # close the total of the reactants and modifiers
        mathXMLString =mathXMLString + """
             </apply>
            """

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
                <cn type="e-notation"> 6 <sep/> 23 </cn>
                <cn type="e-notation"> 5.23 <sep/> -16 </cn>
              </apply>
             </apply>
            """


          if total_conc_stoic == 1:
            mathXMLString =mathXMLString + conc_string
          else:
            mathXMLString =mathXMLString + "<apply><divide/>" +  conc_string + "<cn type=\"integer\">"+ str(total_conc_stoic)+"</cn></apply>"



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
