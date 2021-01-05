This python script generates subnetworks in SBML based on the network of the Mycoplasma genitalium organism (mg_net network).
More information about this reaction network can be found in the paper "A biochemical network modeling of a whole-cell" published by Burke et al.
(Web link: https://www.nature.com/articles/s41598-020-70145-4)
The SBML and graphML files of mg_net network are required in order to be able to run this script.
(Web link for download these files: https://github.com/pauloburke/whole-cell-network).
The graphML file is needed because it includes the information regarding the cellular process of each reaction.
The initial mg_net SBML file does not have any initial values for the model species and any reaction rate formulas for the model reactions.
Therefore, this script also gives initial values to all the species and reaction rate formulas to all the reactions.


usage: generate_network.py [-h] -s <sbml_filename> -g <graphML_filename> [-o <output_filename>] [-p <process_id> [<process_id> ...]] [-r <min> <max>]
                           [-n <mean> <range_value_from_mean>]

required arguments:
  -s <sbml_filename>    Specify the SBML filename of mg_net network
  -g <graphML_filename>
                        Specify the GraphML filename of mg_net network

optional arguments:
  -o <output_filename>  Specify the output filename
  -p <process_id> [<process_id> ...]
                        Specify the types of processes for the new subnetwork (default: all the processes). 
                        Type the process_ids corresponding to the processes separated with space:
                        
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
                          
  -r <min> <max>        Specify the minimum and maximum value for the "discrete uniform" distribution 
                        for the random initialization of species initial count (default from 10 to 1000)
  -n <mean> <range_value_from_mean>
                        Specify the mean and the range value from mean for the normal distribution 
                        for the initialization of k coefficient used in the reaction rate formulas
                        (default mean=5 and range value from mean=2.5)

Examples:

1.The following example generates a subnetwork in SBML with the reactions belonging to the "Replication Initiation" and "Cellular Division" cellular processes.
  
  python3 generate_network.py -s mg_net.sbml -g mg_net.graphml -p 14 15

2.This example generates a network in SBML with all the reactions of mg_net model but with a specific min and max value for 
  the "discrete uniform" distribution for the random initialization of species initial count and save the output to an output file.

  python3 generate_network.py -s mg_net.sbml -g mg_net.graphml -r 100 2000 -o output.txt

3.The following example generates a subnetwork in SBML with only the reactions of the "Replication" cellular process where specific mean 
  and range value from mean for the normal distribution for the initialization of k coefficient is used and the output is saved to an output file. 
  
  python3 generate_network.py -s mg_net.sbml -g mg_net.graphml -p 1 -n 3 1.5 -o output.txt




