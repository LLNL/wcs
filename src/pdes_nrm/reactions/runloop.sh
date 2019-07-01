#!/bin/bash

nlptot=$((16*64*10))

for nn in 1 ; do
    for cpn in 1 3 4 5 8 13 16 ; do
        # Can replace
	#    srun -N$nn -n$(($nn*$cpn)) -t20 ...
	# with
	#    mpirun -np $(($nn*$cpn)) ...
	# for running on workstation or laptop.
	srun -N$nn -n$(($nn*$cpn)) -t20 \
	    ./matrix-cancel-model \
	    --nlp=$(($nlptot/$nn/$cpn)) \
	    --end=$t \
	    --synch=$((($nn*$cpn==1)?1:3)) \
	    --extramem=1000 \
	    --report-interval=0.1 \
	    --batch=2 \
	    --gvt-interval=1024 \
	    --clock-rate=1e6 | tee run.out.$nn.$cpn.t$t
    done
done

# These lines print each output file name (for each
# node/MPI rank configuration). It also does a diff
# with data from serial run for correctness checking.
# If there is any output from these diff's something
# is broken.
for f in run.out.* ; do
    echo "### $f ###"
    diff <(grep -A 5 '^  A =' run.out.1.1.t100) <(grep -A 5 '^  A =' $f)
done
