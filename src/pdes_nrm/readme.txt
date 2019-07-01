This repository contains Tomas' version of ROSS, used in e.g.
his KMC code, and a model for testing. This model allows checking
correctness in a more discriminating way than PHOLD.

Subdirectories:

  matrix: A model used for debugging. It is kind of like PHOLD, but
          also enables proving that messages were received in order,
	  using a finite ring matrix algebra.

  ross:   Tomas' version of ROSS

  matrix-cancel: A model used for debugging and illustration of event
                 cancellation (retraction). Also exercises commit
		 function.

Building the code:

1. Build ROSS:
   cd ross/src
   # Edit Makefile to set copmiler name (one of the first lines)
   # currectly the compiler is mpicxx (MPI C++ compiler necessary).
   # Make and copy library to ../lib directory:
   make
   make install
   cd ../..

2. Build matrix model:
   cd matrix
   # Set compiler also in compile.matrix-model
   # Build:
   . compile.matrix-model
   cd ..
   
3. Run test:
   cd matrix
   # Test run script usus mpirun to launch a parallel job.
   # this is a small test and runs as is on many login nodes,
   # workstations, and laptops. For other machines, getting
   # a node allocation might be required first, and perhaps
   # replacing mpirun with e.g. srun.
   ./run.matrix-model

   # The output of the test run can be compared for correctness
   # with sample output at the end of the run.matrix-model file.
   # See that file for more details.
   cd ..

4. Build matrix-cancel model:
   cd matrix-cancel
   . compile.matrix-cancel
   cd ..

5. Run matrix-cancel test:
   cd matrix-cancel
   ./runloop.sh

   Test takes maybe a couple of minutes. Lots of text data is eventually
   dumped to the screen. This data is also saved in files with names
   run.out.x.y where x is number of nodes (1 in runloop.sh at this time),
   and y is number of MPI tasks per node.

   At the end of the test, the name of each run.out.* file is printed
   in between #'s. After each name is printed is a diff of relevant
   data betwen this output file and a sequential run. If the runs
   (code) are correct, there should be no output from diff, and just
   the file names are listed.


4. Enjoy:-)
