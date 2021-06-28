## Python script for generating *Mycomplasma genitalium* networks
This python script generates subnetworks in SBML based on the network of the Mycoplasma genitalium organism (mg_net network).
More information about this reaction network can be found in Burke, P.E.P., Campos, C.B.d.L., Costa, L.d.F. et al. [A biochemical network modeling of a whole-cell.](https://www.nature.com/articles/s41598-020-70145-4) Sci Rep 10, 13303 (2020).
The SBML and graphML files of mg_net network are required to run this script
([Download these files](https://github.com/pauloburke/whole-cell-network)).
The GraphML file contains the information about which cellular process each reaction belongs to. We use this information to select reactions in the SBML file that belongs to the cellular processes chosen by a user and generate an SBML file of the sub-network that only contains those reactions.
The network file does not include any information regarding kinetic reaction rates or initial concentrations of the species, both of which are required to simulate the model.
Therefore, this script also assigns an initial population to each species and a reaction rate formula to each reaction.
### Required packages:
+ **argparse**
+ **xml**
+ **libsbml**

### Usage:
`generate_network.py [-h] -s <sbml_filename> -g <graphML_filename> [-o <output_filename>] [-p <process_id> [<process_id> ...]] [-u <min> <max>] [-n <min> <max>]`

**Required arguments and inputs**:

  `-s <sbml_filename>`:    Specify the SBML filename of mg_net network ([mg_net.sbml](https://github.com/pauloburke/whole-cell-network))

  `-g <graphML_filename>`:
                        Specify the GraphML filename of mg_net network ([mg_net.graphml](https://github.com/pauloburke/whole-cell-network))

**Optional arguments**:

  `-o <output_filename>`: Specify the output filename

  `-p <process_id> [<process_id> ...]`:
                        Specify the cellular processes to be included in the new subnetwork (default: all the processes) by listing the process_ids separated by a space. 
                        
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
                          
  `-u <min> <max>`:        Specify the range of the "discrete uniform" distribution
used in the random initialization of the initial copy-number of each species (default from 30000 to 500000)

  `-n <min> <max>`:
                        Specify the range of the normal distribution 
                        for randomly assigning reaction rate coefficients used in rate formulas
                        (default from 2.5 to 7.5). The mean of the normal distribution is calculated as (max-min)/2. The standard deviation is approximated by (max-mean)/3 to cover the 99.7% confidence interval within min and max based on the three-sigma rule of thumb. Finally, we replace any out of the range value by redrawing to make sure all of random values fall within the range. k is the same for all the reactions in a cellular process.

  `-i <input_count_species_filename>`:
                        Specify the SBML filename for the initial count of the species for the metabolism cellular 
                        process species contained in the specific network.

  `-x <approach_id>`:
                        Specify which approach you want to generate (1 for the first approach or 2 for the second approach
                        or 3 for the third approach or 23 for the combined approach 2 and 3 )

  `-m <max_sum_stoichiometry>`:
                        Specify the maximum sum of stoichiometries as the threshold between the first and third approach. 
                        Value should be between 3 and 10 (default value 7).
                        

### Examples:

1. The following example generates a subnetwork in SBML with the reactions of the cellular processes, "Replication Initiation" and "Cellular Division".
  ```
  python3 generate_network.py -s mg_net.sbml -g mg_net.graphml -p 14 15
  ```
2. This example generates a network in SBML including all the reactions of mg_net model while randomly initializing the species count using a uniform distribution between min and max. The output is saved to an output file.
  ```
  python3 generate_network.py -s mg_net.sbml -g mg_net.graphml -u 100 2000 -o output.txt
  ```
3. The following example generates a subnetwork in SBML including only the reactions of the "Replication" cellular process while randomly assigning the reaction rate coefficient using a normal distribution between min and max. The output is saved to an output file. 
  ```
  python3 generate_network.py -s mg_net.sbml -g mg_net.graphml -p 1 -n 1.5 4.5 -o output.txt
  ```




