#!/bin/bash
# This script takes a SSA trajectory result file as an input, and prints out
# the trajectory of a chosen subset of species (or reactions) in the result.
# The species (or reactions) to print out are selected by the other script
# `./get_max_species_list.awk`.
# It finds the top-N species (or reactions) in terms of the initial copy number,
# the final copy number, and the difference between the initial and the final
# copy numbers.

if [ $# -ne 1 ] ; then
  echo usage: $0 trajectory_file
  exit 0
fi

trj_file=$1

show_reactions=0

num_max_initial=1
num_max_diff=0
num_max_dyn=10
num_max_final=1

if [ "${show_reactions}" == "1" ] ; then
  num_max_initial=0
  if [ ${num_max_final} -lt ${num_max_dyn} ] ; then
    num_max_final=${num_max_dyn}
  fi
  num_max_dyn=0
  if [ ${num_max_final} -lt ${num_max_diff} ] ; then
    num_max_final=${num_max_diff}
  fi
  num_max_diff=0
fi

list=`./get_max_species_list.awk \
  -v num_max_initial=${num_max_initial} \
  -v num_max_diff=${num_max_diff} \
  -v num_max_dyn=${num_max_dyn} \
  -v num_max_final=${num_max_final} \
  -v show_reactions=${show_reactions} \
  ${trj_file} | awk \
  '{
    if (($1=="Unique") && ($2 == "indices:")) {
      printf("%u", $3);
      for(i=4; i <= NF; i++) {
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
