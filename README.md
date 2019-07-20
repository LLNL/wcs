# Whole Cell Simulator
 A whole cell model (WCM) is a comprehensive multi-scale computational model
 representing all the known biochemical processes in a cell. It relies on a
 variety of intracellular pathway models and omics data.

 This computational framework enables seamless integration of diverse simulation
 methods used in WCM such as stochastic simulation algorithm (SSA), ordinary
 differential equations (ODE), flux balance analysis (FBA), and logic-based
 approaches.
 These methods run simultaneously, not only for whole pathways but also for
 subsets of reactions. Furthermore, we allow dynamically switching between
 methods when beneficial.

## Current Requirements:
 + c++ compiler that supports c++17.
   e.g., clang++ 5.0 or later, g++ 7.1 or later, and icpc 19.0.1 or later
 + GNU Boost library
   Currently, there are three versions of boost available on quartz.
   The one under the default system path, another under
   "/usr/tce/packages/boost/boost-1.69.0-mvapich2-2.2-gcc-4.9.3", and finally
   the one under "/usr/tce/packages/boost/boost-1.66.0-mvapich2-2.2-gcc-6.1.0"

   The first one does not work well with the compiler choices above.
   The second one works well with the clang choice, and the third with gcc.
   To override building with a manual choice of boost, set the environment
   variable `BOOST_ROOT` or pass `-DBOOST_ROOT=<path-to-the-chosen-boost-build>`
   to cmake. To run the executable, you might need to add `${BOOST_ROOT}/lib`
   to the `LD_LIBRARY_PATH`
 + cmake 3.12 or later

## Future requirements:
 + Charm++ and Charades (ROSS over Charm++)
 + Sundial CVODE

## Optional requirements:
 + ExprTK
 + libSBML
