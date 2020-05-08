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

## Getting started:
 ```
 git clone https://github.com/llnl/wcs.git
 mkdir build; cd build
 cmake -DBOOST_ROOT:PATH=<PathToYourBoostDev> -DCMAKE_INSTALL_PREFIX:PATH=<YourInstallPath> ../wcs
 make
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

 We rely on either of the two approaches to ingest problem inputs. The one
 based on libSBML, currently under development. The other is based on ExprTk
 to parse the formula of the reaction rate annotated to each reaction
 vertex in the GraphML-formatted description of the reaction network.
 To enable one of the optional components, either set the environment
 variable `WCS_WITH_SBML` or `WCS_WITH_EXPRTK`, or invoke cmake with
 the option `-DWCS_WITH_SBML:BOOL=ON` or `-DWCS_WITH_EXPRTK:BOOL=ON`
 respectively for libSBML or ExprTk.
 In case that both are set, only libSBML is enabled and ExprTk is ignored.
 To use a pre-installed libSBML, either set the environment variable
 `SBML_ROOT` or invoke cmake with the option `-DSBML_ROOT=<path-to-libsbml>`.
 For ExprTk, no pre-installation is required as it will automatically download
 it during the build if the presence is not detected. To use a pre-installed
 copy, use the variable `EXPRTK_ROOT`.

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
