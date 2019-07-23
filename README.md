Whole Cell Simulator

Current Requirements:
 - c++ compiler that supports c++17.
   e.g., clang++ 5.0 or later, g++ 7.1 or later, and icpc 19.0.1 or later
 - GNU Boost library
   Currently, there are three versions of boost available on quartz.
   The one under the default system path, another under
   "/usr/tce/packages/boost/boost-1.69.0-mvapich2-2.2-gcc-4.9.3", and finally
   the one under "/usr/tce/packages/boost/boost-1.66.0-mvapich2-2.2-gcc-6.1.0"

   The first one does not work well with the compiler choices above.
   The second one works well with the aforementioned clang version, and the
   third with that of gcc.
   To override building with a manual choice of boost, set the environment
   variable "BOOST\_ROOT" or pass "-DBOOST\_ROOT=<path-to-the-chosen-boost-build>"
   to cmake. To run the executable, you might need to add "${BOOST\_ROOT}/lib"
   to the "LD\_LIBRARY\_PATH"

Future requirements:
 - Charm++ and Charades (ROSS over Charm++)
 - Sundial CVODE

Optional requirements:
 - ExprTK
 - libSBML
