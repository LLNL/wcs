#!/bin/awk -f
# This script takes an SSA trajectory result file, and finds the top-N
# species in terms of the initial copy number, the final copy number,
# the difference between the initial and the final copy numbers, the
# sum of absolute differences between subsequent data points.
#
# > ./get_max_species_list.awk [options] trajectoy_file
#
# By default, N is set to 10. To change this, use the following options
#   -v num_max_initial=N1
#   -v num_max_diff=N2
#   -v num_max_dyn=N2
#   -v num_max_final=N3

# The format of the trajectory result file is as follows.
# First line shows the number of species, reactions and events.
# Second line shows labels of columns. The first column is the time of event
# with label "Time". The rest is species names. Each of the rest lines
# corresponds to an event in the chronological order.

function abs(v) {return int(v < 0 ? -v : v);}

(NR==1) {
  if (length(num_max_initial) == 0) {
    num_max_initial = 10;
    for (j=1; j <= num_intial_max; j++) {
      max_init_val[j] = -1;
    }
  }
  if (length(num_max_diff) == 0) {
    num_max_diff = 10;
    for (j=1; j <= num_max_diff; j++) {
      max_diff_val[j] = -1;
    }
  }
  if (length(num_max_dyn) == 0) {
    num_max_dyn = 10;
    for (j=1; j <= num_max_dyn; j++) {
      max_dyn_val[j] = -1;
      prev_val[j] = 0;
      sum_diffs[j] = 0;
    }
  }
  if (length(num_max_final) == 0) {
    num_max_final = 10;
    for (j=1; j <= num_max_final; j++) {
      max_final_val[j] = -1;
    }
  }
}

(FNR>1) {
  if (FNR == 2) {
    for (i=2; i <= NF; i++) {
      species[i] = $i # Store species names
    }
    next;
  } else if (FNR == 3) {
    for (i=2; i <= NF; i++) {
      init_val[i] = $i; # Store initial population
      # Find the top `num_max_initial` species in terms of the initial population
      val = $i;
      idx = i;
      for (j=1; j <= num_max_initial; j++) {
        if (max_init_val[j] < val) {
          tmp_val = max_init_val[j];
          max_init_val[j] = val;
          val = tmp_val;
          tmp_idx = max_init_idx[j];
          max_init_idx[j] = idx;
          idx = tmp_idx;
        }
      }
      prev_val[i] = $i;
    }
  } else {
    for (i=2; i <= NF; i++) {
      sum_diffs[i] += abs($i - prev_val[i]);
    }
  }
}

END {
  # Find the species with the maximum difference between the initial
  # and the final population
  for (i=2; i <= NF; i++) {
    val = abs($i - init_val[i]);
    idx = i;
    for (j=1; j <= num_max_diff; j++) {
      if (max_diff_val[j] < val) {
        tmp_val = max_diff_val[j];
        max_diff_val[j] = val;
        val = tmp_val;
        tmp_idx = max_diff_idx[j];
        max_diff_idx[j] = idx;
        idx = tmp_idx;
      }
    }
  }
  # Find the top `num_max_final` species in terms of the final population
  for (i=2; i <= NF; i++) {
    val = $i;
    idx = i;
    for (j=1; j <= num_max_final; j++) {
      if (max_final_val[j] < val) {
        tmp_val = max_final_val[j];
        max_final_val[j] = val;
        val = tmp_val;
        tmp_idx = max_final_idx[j];
        max_final_idx[j] = idx;
        idx = tmp_idx;
      }
    }
  }

  k = 0;
  printf("initial max :");
  for (i=1; i <= num_max_initial; i++) {
    printf("  %s %d", species[max_init_idx[i]], max_init_val[i]);
    #printf("  %s", species[max_init_idx[i]]);
    if (exists[max_init_idx[i]] == 0) {
      exists[max_init_idx[i]] = 1;
      uniq[++k] = max_init_idx[i];
    }
  }
  printf("\n");

  printf("max diff    :");
  for (i=1; i <= num_max_diff; i++) {
    printf("  %s %d", species[max_diff_idx[i]], max_diff_val[i]);
    #printf("  %s", species[max_diff_idx[i]]);
    if (exists[max_diff_idx[i]] == 0) {
      exists[max_diff_idx[i]] = 1;
      uniq[++k] = max_diff_idx[i];
    }
  }
  printf("\n");

  printf("max dynamics:");
  for (i=2; i <= NF; i++) {
    # Find the top `num_max_dyn` species in terms of the sum of the
    # difference between subsequent data points
    val = sum_diffs[i];
    idx = i;
    for (j=1; j <= num_max_dyn; j++) {
      if (max_dyn_val[j] < val) {
        tmp_val = max_dyn_val[j];
        max_dyn_val[j] = val;
        val = tmp_val;
        tmp_idx = max_dyn_idx[j];
        max_dyn_idx[j] = idx;
        idx = tmp_idx;
      }
    }
  }
  for (i=1; i <= num_max_dyn; i++) {
    printf("  %s %d", species[max_dyn_idx[i]], max_dyn_val[i]);
    #printf("  %s", species[max_dyn_idx[i]]);
    if (exists[max_dyn_idx[i]] == 0) {
      exists[max_dyn_idx[i]] = 1;
      uniq[++k] = max_dyn_idx[i];
    }
  }
  printf("\n");

  printf("max final   :");
  for (i=1; i <= num_max_final; i++) {
    printf("  %s %d", species[max_final_idx[i]], max_final_val[i]);
    #printf("  %s", species[max_final_idx[i]]);
    if (exists[max_final_idx[i]] == 0) {
      exists[max_final_idx[i]] = 1;
      uniq[++k] = max_final_idx[i];
    }
  }
  printf("\n");

  asort(uniq, unique);

  printf("The number of unique species: %d\n", k);
  printf("Unique species indices:");
  for (i=1; i <= k; i++) {
    printf(" %d", unique[i]);
  }
  printf("\n");

  printf("Unique species names:");
  for (i=1; i <= k; i++) {
    printf(" %s", species[unique[i]]);
  }
  printf("\n");
}
