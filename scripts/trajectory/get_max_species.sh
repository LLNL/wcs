#!/bin/bash
# This script takes a SSA trajectory result file as an input, and prints out
# the trajectory of a subset of species in the result. The species to print
# out are selected by the other script `./get_max_species_list.awk`.
# It finds the top-N species in terms of the initial copy number, the final
# copy number, and the difference between the initial and the final copy
# numbers.

if [ $# -ne 1 ] ; then
  echo usage: $0 trajectory_file
  exit 0
fi

trj_file=$1

num_max_initial=10
num_max_diff=10
num_max_dyn=10
num_max_final=10

list=`./get_max_species_list.awk \
  -v num_max_initial=${num_max_initial} \
  -v num_max_diff=${num_max_diff} \
  -v num_max_dyn=${num_max_dyn} \
  -v num_max_final=${num_max_final} \
  ${trj_file} | awk \
  '{
    if (($1=="Unique") && ($2 == "species") && ($3 == "indices:")) {
      printf("%u", $4);
      for(i=5; i <= NF; i++) {
        printf(" %u", $i);
      }
      printf("\n");
    }
  }'`


awk -v list="$list" \
'(NR == 1) {
  n = split(list, indices," ");
}
(NR == 2) {
  printf("Time");
  for(i = 1 ; i <= n ; i++) {
    #printf("%s\n", indices[i]);
    printf(" %s", $(indices[i]));
  }
  printf("\n");
}
(NR > 2) {
  printf("%s", $1);
  for(i = 1 ; i <= n ; i++) {
    printf(" %s", $(indices[i]));
  }
  printf("\n");
}' ${trj_file}
