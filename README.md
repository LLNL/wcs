# Whole Cell Simulator
 A whole cell model (WCM)  is a comprehensive multi-scale computational
 model representing all the known biochemical processes in a cell. It relies
 on a variety of intracellular pathway models and omics data.

 This computational framework  enables  seamless integration of diverse
 simulation methods used in WCM such as stochastic simulation algorithm
 (SSA), ordinary differential equations (ODE), flux balance analysis (FBA),
 and logic-based approaches.
 These methods run simultaneously, not only for whole pathways but
 also for subsets of reactions. Furthermore, we aim to enable switching
 dynamically between the methods when beneficial.

## Current Requirements:
 + **Platforms targeted**: Linux-based systems
 + **c++ compiler that supports c++17**
   e.g., clang++ 5.0 or later, g++ 7.1 or later, and icpc 19.0.1 or later

 + **GNU Boost library**
   The particular boost libraries used are `graph`, `filesystem`, `regex` and
   `system`.
   To build with other pre-existing installation of boost, set the environment
   variable `BOOST_ROOT` or pass `-DBOOST_ROOT=<path-to-boost>`
   to cmake. An example path to boost on LC is
   `/usr/tce/packages/boost/boost-1.69.0-mvapich2-2.3-gcc-8.1.0`.
   To run the executable, add `${BOOST_ROOT}/lib` to the `LD_LIBRARY_PATH` as
   needed

 + **Guide specific to using clang on Livermore Computing (LC) platforms**
   Currently, to use clang on Livermore Computing (LC) platforms with
   the libraries compiled with gcc, make sure the c++ standard library
   is compatible. On LC, clang by default is paired with c++ standard
   library from gcc/4.9.3. To avoid incompatibility issue,

   0) Use the user libraries compiled with the system default compiler.
      On many of LC platforms, it is currently gcc/4.9.3. On others,
      it depends on how clang is configure there.
   1) Make clang use the c++ standard library from the same version of gcc
      as that used for building user libraries to link.
      e.g., clang++ --gcc-toolchain=/usr/tce/packages/gcc/gcc-8.1.0/ ...
   2) Use clang's c++ standard library. Recompile user libraries in the
      same way as needed.
      i.e., clang++ --stdlib=libc++ ...

   Choose either `USE_GCC_LIBCXX` for option 1 or `USE_CLANG_LIBCXX` for
   option 2 if needed. Usually, option 1 works best. Option 0 does not work
   well especially with Catch2 due to the incomplete support for C++17.
   If neither is chosen, the build relies on the system default, which is,
   on LC, with `USE_GCC_LIBCXX` on and `GCC_TOOLCHAIN_VER` set to "4.9.3".
   If both are on, `USE_GCC_LIBCXX` is turned off. When `USE_GCC_LIBCXX`
   is on, `GCC_TOOLCHAIN_VER` can be set accordingly (e.g., "8.1.0").

 + **Guide specific to using the intel compiler on Livermore Computing (LC) platforms**
   For gcc Interoperability, use `-DGCC_PATH=<path-to-gcc>`. Otherwise,
   icc picks up the gcc in default path, which may not support c++17.

 + **cmake 3.12 or later**
   This requirement mostly comes from the compatibility between the cmake
   module `find_package()`, and the version of boost used. An older version
   might still work with some warnings.

 + [**Cereal**](https://uscilab.github.io/cereal)
   We rely on Cereal serialization library to enable state packing and
   unpacking of which needs arises under various circumstances including
   messaging, rollback, migration, and checkpointing. Cereal is a c++
   header-only library. No pre-installation is required as it is
   automatically downloaded and made available.

 + [**iheap**](https://github.com/yangle/iheap.git)
   We rely on iheap, the indexed-heap library, to implement the next reaction
   method. This is a header-only library consists of a single file.
   No pre-installation is required as it is automatically downloaded and made
   available.

 + [**Protocol Buffers**](https://developers.google.com/protocol-buffers)
   We use the google protocol buffers library for parsing the configuration file
   of simulation, which is written by users in the [**protocol buffers language**](https://developers.google.com/protocol-buffers/docs/proto3).
   This is a required package. A user can indicate the location of a
   pre-installed copy via `-DPROTOBUF_ROOT=<path>`. Without it, building WCS
   consists of two stages. In the first stage, the source of protocol buffer will
   be downloaded. Then, the library as well as the protoc compiler will be built
   and installed under where the rest of WCS project will be.
   In the second stage, the WCS project will be built using the protocol buffer
   installed in the first stage. Both stages require the same set of options for
   the cmake command.
   In case of cross-compiling, the path to the protoc compiler and the path to
   the library built for the target platform can be explicitly specified via
   `-DProtobuf_PROTOC_EXECUTABLE=<installation-for-host/bin/protoc>`
   and `-DPROTOBUF_DIR=<installation-for-target>` respectively.

## Getting started:
 ```
 git clone https://github.com/llnl/wcs.git
 mkdir build; cd build
 # The first invocation of cmake will setup to build protocol buffer
 cmake -DBOOST_ROOT:PATH=<PathToYourBoostDev> \
       -DCMAKE_INSTALL_PREFIX:PATH=<YourInstallPath> \
       -DWCS_WITH_SBML:BOOL=ON \
       -DSBML_ROOT:PATH=<path-to-libsbml> \
       ../wcs
 make -j 4
 # The second invocation of cmake will setup to build the WCS project
 cmake -DBOOST_ROOT:PATH=<PathToYourBoostDev> \
       -DCMAKE_INSTALL_PREFIX:PATH=<YourInstallPath> \
       -DWCS_WITH_SBML:BOOL=ON \
       -DSBML_ROOT:PATH=<path-to-libsbml> \
       ../wcs
 make -j 4
 make install
 ```

## Future requirements:
 + **Charm++ and Charades (ROSS over Charm++)**
 + [**Sundial CVODE**](https://github.com/LLNL/sundials.git)
   For linking with CVODES of Sundials package, use `-DWCS_WITH_SUNDIALS:BOOL=ON`
   and `-DSUNDIALS_ROOT:FILEPATH=<path-to-sundials>`.
   Make sure that Sundials is built with the cmake option `-DBUILD_CVODE:BOOL=ON`.

## Optional requirements:
 + [**libSBML C++ API**](http://sbml.org/Software/libSBML)
 + [**ExprTK**](https://github.com/ArashPartow/exprtk)

	While these are optional, either of the two frameworks or both must be used.

	WCS can rely on SBML for input model ingestion. An input model contains MathML-
 based description of reaction formula. WCS relies on libSBML to parse the
 formula, then generates a c++ code that includes callable functions for the
 formula, which is then compiled and linked on-line.

	Alternatively, WCS can use ExprTk to parse the formula of the reaction rate
 annotated to each reaction in a GraphML-formatted input. ExprTk dynamically
 generates a callable function for each reaction formula described in the input
 model. For GraphML-based inputs, users must express reaction rate formula in
 a syntax compatible with ExprTk.

	When both SBML and ExprTk are enabled, WCS relies on libSBML to parse the
 SBML-based input models. Then, WCS internally generates formula strings that
 are compatible with ExprTk to take advantage of it for function generation.
 With SBML enabled but without ExprTk, on-line code-generation scheme is used.
 With ExprTk but without SBML, GraphML-based input ingestion scheme is used.

	To enable SBML-based ingestion, set the environment variable `WCS_WITH_SBML`
 or invoke cmake command with the option `-DWCS_WITH_SBML:BOOL=ON`.
 To enable ExprTk-base function generation for reaction formula, use the
 environment variable `WCS_WITH_EXPRTK`, or use the cmake option
 `-DWCS_WITH_EXPRTK:BOOL=ON`

	To use a pre-installed libSBML, either set the environment variable
 `SBML_ROOT` or invoke cmake with the option `-DSBML_ROOT=<path-to-libsbml>`.
 For ExprTk, no pre-installation is required as it will automatically download
 it during the build if the presence is not detected. To use a pre-installed
 copy, use the variable `EXPRTK_ROOT`.

 + [**Metis**](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
 WCS rely on Metis graph partitioning library to perform initial partitioning
 of reaction network for parallel execution. The partitioning outcome will
 further be refined by custom algorithm. To link with Metis, use the cmake
 options `-DWCS_WITH_METIS=ON` and `-DMETIS_ROOT=<path-to-metis-install>`.

## Taking advanage of OpenMP:
 Build WCS with the cmake option `-DWCS_WITH_OPENMP:BOOL=ON`.
 Then, for accelerating execution performance, use the OpenMP environment
 variables to control the parallelism and the processor affinity, such as
 `OMP_NUM_THREADS`, `OMP_PROC_BIND`, and `OMP_PLACES`.
 Currently, only the next reaction method is parallelized using OpenMP.
 There are multiple regions where OpenMP parallelization is applied.
 Use the cmake option `-DWCS_OMP_MODE=[1-5]` to enable parallelization of a
 particular combination of regions. The mode 1, 2, and 3 enables fine-grained
 parallelizaion. Mode 1 helps when there exists a large number of
 interdependent reactions. In other words, firing a reaction leads to
 updating a large number of reactions dependent on the species updated by
 the reaction fired. Mode 2 is similar to mode 1 and helps additional cases
 with a large number of reactants per reactions. Mode 3 further includes
 additional cases of a large number of products per reaction.
 The mode 4, and 5 implements the parallelization using graph partitioning.
 These two modes are likely to be suitable for a broader range of problems,
 especially when the size of network is not small. Metis is required for these.

## Unit testing:
 + [**Catch2**](https://github.com/catchorg/Catch2)
 We rely on Catch2 for unit testing. To enable testing for development, set
 the environment variable `WCS_WITH_UNIT_TESTING` or invoke cmake with the
 option `-DWCS_WITH_UNIT_TESTING=ON`. No pre-installation of Catch2 is required
 as it is automatically downloaded and made available. To use a pre-installed
 copy, use the variable `CATCH2_ROOT`.

## Authors:
  Many thanks go to WCS's [contributors](https://github.com/llnl/wcs/graphs/contributors).

 * Jae-Seung Yeom
 * Konstantia Georgouli
 * Giorgis Georgakoudis
 * Robert Blake
 * Ali Navid

## Release:
 Whole cell simulator is distributed under the terms of the MIT license.
 All new contributions must be made under this license.
 See [LICENSE](https://github.com/llnl/wcs/blob/master/LICENSE) and [NOTICE](https://github.com/llnl/wcs/blob/master/NOTICE) for details.

 + `SPDX-License-Identifier: MIT`
 + `LLNL-CODE-809742`

## Contributing:
 Please submit any bugfixes or feature improvements as [pull requests](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork).
