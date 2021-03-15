##### `get_max_species_list.awk`

This script takes an SSA trajectory result file, and finds the top-`N` species
or reactions in terms of the initial copy number, the final copy number, the
the difference between the initial and the final copy numbers, the sum of
absolute differences between subsequent data points.
 ```
 > ./get_max_species_list.awk [options] trajectoy_file
 ```
 By default, `N` is set to 10. To change this, use the following options:
 ```
   -v num_max_initial=N1
   -v num_max_diff=N2
   -v num_max_dyn=N3
   -v num_max_final=N4
 ```
 The result returned is the combined set of elements that satisfy each
 criteria. Each element in the result is represented by either the name
 or the column id.

 By default, species statistics are shown. Use the following option to get
 reaction statistics
 ```
   -v show_reactions=1
 ```

 As the number of firing of a reaction monotonically increases, all the
 selection criteria other than `num_max_initial` is the same for reactions.
 Also note that the initial count of reactions fired are zero for every
 reaction, and thus `num_max_initial` is not a meaningful criterion for
 selecting reactions.

 The format of the trajectory result file is as follows. The first line shows
 the number of species, reactions and events. The second line shows labels of
 columns. The first column is the time of event with label "Time". The rest is
 species names followed by reaction names. Each of the rest lines corresponds
 to an event in the chronological order.

Example:
```
$ ./get_max_species_list.awk -v num_max_initial=2 -v num_max_final=2 -v num_max_diff=0 trajectory.txt
initial max :  c_GTP 2000  c_TranslationComplex_MG_316_MONOMER 2000
max diff    :
max dynamics:  c_GDP 305106  c_PI 305106  c_GTP 29619  c_MG501 24037  c_MG514 22135  c_MG472 19924  c_MG511 19179  c_MG490 18506  c_MG513 17450  c_MG500 16472
max final   :  c_GDP 77061  c_PI 76930
The number of unique species: 11
Unique indices: 2 3 486 506 520 522 540 544 546 773 1186
Unique names: c_GDP c_GTP c_MG472 c_MG490 c_MG500 c_MG501 c_MG511 c_MG513 c_MG514 c_PI c_TranslationComplex_MG_316_MONOMER
```


#### `./get_max_species.sh`

This script takes an SSA trajectory result file as an input. Then, it prints out
the trajectory of a subset of species (or reactions). It relies on the other
script `./get_max_species_list.awk` to select the species (or reactions) to print
out.  For adjusting the number of species (or reactions) to select, change the
variables defined in the script.
