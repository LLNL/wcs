cmake_minimum_required(VERSION 3.12)

project(WCS CXX)

# Prevent in-source builds
if (PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(FATAL_ERROR
    "In-source build attempted; please clean the CMake cache and then "
    "switch to an out-of-source build, e.g.,\n"
    "rm -rf CMakeCache.txt CMakeFiles/\nmkdir build && "
    "cd build && cmake <options> ..\n")
endif ()

# Add CMake modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
include(WCSCMakeUtilities)

#FIXME
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if (CMAKE_BUILD_TYPE MATCHES Debug)
  set(WCS_DEBUG TRUE)
endif ()

if (NOT DEFINED BUILD_SHARED_LIBS)
  set(BUILD_SHARED_LIBS ON)
endif ()

# Convenience defines
string(TOUPPER "${PROJECT_NAME}" UPPER_PROJECT_NAME)
string(TOLOWER "${PROJECT_NAME}" LOWER_PROJECT_NAME)


################################################################
# Suppress warning
# "Manually-specified variables were not used by the project"
################################################################

set(ignoreUnused ${WCS_OMP_MODE}${WCS_WITH_VTUN}${VTUNE_DIR})
set(ignoreUnused ${WCS_WITH_HPCTOOLKIT}${HPCTOOLKIT_DIR})
set(ignoreUnused ${WCS_WITH_ROSS}${ROSS_ROOT}${ROSS_DIR})

################################################################
# Version setup
################################################################

set(WCS_VERSION_MAJOR 0)
set(WCS_VERSION_MINOR 1)
set(WCS_VERSION_PATCH 0)

# Build with at least C++17 standard; allow newer standards.
if (NOT CMAKE_CXX_STANDARD OR CMAKE_CXX_STANDARD EQUAL 98
    OR CMAKE_CXX_STANDARD EQUAL 11)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
  set(CMAKE_CXX_EXTENSIONS ON)
endif ()

# This will help define some targets later
if (CMAKE_VERSION VERSION_LESS 3.9)
  set(WCS_PUBLIC_LINK_FLAG)
else ()
  set(WCS_PUBLIC_LINK_FLAG "PUBLIC")
endif ()

set(WCS_VERSION "${WCS_VERSION_MAJOR}.${WCS_VERSION_MINOR}.${WCS_VERSION_PATCH}")

# Check to see if we are in a git repo
find_program(__GIT_EXECUTABLE git)
mark_as_advanced(__GIT_EXECUTABLE)
if (__GIT_EXECUTABLE)

  execute_process(
    COMMAND ${__GIT_EXECUTABLE} rev-parse --is-inside-work-tree
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    OUTPUT_VARIABLE __BUILDING_FROM_GIT_SOURCES
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  if (__BUILDING_FROM_GIT_SOURCES)
    # Get the git version so that we can embed it into the executable
    execute_process(
      COMMAND ${__GIT_EXECUTABLE} rev-parse --show-toplevel
      WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
      OUTPUT_VARIABLE __GIT_TOPLEVEL_DIR
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(
      COMMAND ${__GIT_EXECUTABLE} rev-parse --git-dir
      WORKING_DIRECTORY "${__GIT_TOPLEVEL_DIR}"
      OUTPUT_VARIABLE __GIT_GIT_DIR
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(
      COMMAND ${__GIT_EXECUTABLE} --git-dir "${__GIT_GIT_DIR}" describe
      --abbrev=7 --always --dirty --tags
      WORKING_DIRECTORY "${__GIT_TOPLEVEL_DIR}"
      OUTPUT_VARIABLE __GIT_DESCRIBE_VERSION
      OUTPUT_STRIP_TRAILING_WHITESPACE)

    set(WCS_GIT_VERSION "${__GIT_DESCRIBE_VERSION}"
      CACHE STRING "WCS's version string as told by git.")
  endif (__BUILDING_FROM_GIT_SOURCES)
endif (__GIT_EXECUTABLE)

if (CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
  set(WCS_GNU_LINUX TRUE)
endif ()

################################################################
# Options
################################################################

option(WCS_64BIT_CNT "Enable 64-bit species counter. The default is 32-bit." OFF)

# Sundials may become requirement later
option(WCS_WITH_SUNDIALS "Enable SUNDIALS library" OFF)

option(WCS_WITH_SBML "Enable SBML library" ON)

option(WCS_WITH_EXPRTK "Include ExprTk" OFF)

# Cereal is going to be required.
option(WCS_WITH_CEREAL "Include Cereal" ON)

option(WCS_WITH_VTUNE
  "Statically link the Intel VTune profiling library" OFF)

option(WCS_WITH_HPCTOOLKIT
  "Dynamically link the HPCToolkit profiling library" OFF)

option(WCS_WITH_ROSS "Use ROSS" OFF)

option(WCS_WITH_OPENMP "Enable OpenMP" OFF)

option(WCS_WITH_METIS "Enable METIS library" OFF)

option(WCS_WITH_UNIT_TESTING
  "Enable the unit testing framework (requires Catch2)" OFF)

option(WCS_DEBUG_PRINT_SUBTARGETS
  "Turn on debugging output of internal target properties." OFF)
mark_as_advanced(WCS_DEBUG_PRINT_SUBTARGETS)

# This option is off by default because non-developers should not use
# this option under normal circumstances.
option(WCS_WARNINGS_AS_ERRORS
  "Build with warnings promoted to errors." OFF)
mark_as_advanced(WCS_WARNINGS_AS_ERRORS)

if (NOT WCS_WITH_SBML AND NOT WCS_WITH_EXPRTK)
  message(FATAL_ERROR
          "Requires either or both of SBML and ExprTk.\n")
endif ()

#
# The BGL adjacency list container selector. By default, vecS is used
# for both out-edge lists and vertex lists.
#
#set(WCS_VERTEX_LIST_TYPE "::boost::listS"
#  CACHE STRING "The BGL vertex list type to use in WCS")

######################### For clang ###############################
#set(WCS_OUT_EDGE_LIST_TYPE "::boost::vecS"
#  CACHE STRING "The BGL out edge list type to use in WCS")

#  Setup specific to using clang on Livermore Computing (LC) platforms
#
#  Currently, to use clang on Livermore Computing (LC) platforms with
#  the libraries compiled with gcc (e.g. boost), make sure the c++ standard
#  library is compatible. On LC, clang by default is paired with c++ standard
#  library from gcc/4.9.3. To avoid incompatibility issue,
#
#  0) Use the user libraries compiled with the system default compiler.
#     On many of LC platforms, it is currently gcc/4.9.3. On others,
#     it depends on how clang is configure there.
#  1) Make clang use the c++ standard library from the same version of gcc
#     as that used for building user libraries to link.
#     e.g., clang++ --gcc-toolchain=/usr/tce/packages/gcc/gcc-8.1.0/ ...
#  2) Use clang's c++ standard library. Recompile user libraries in the
#     same way as needed.
#     i.e., clang++ --stdlib=libc++ ...
#
#  Choose either `USE_GCC_LIBCXX` for option 1 or `USE_CLANG_LIBCXX` for
#  option 2 if needed. Usually, option 1 works best. Option 0 does not work
#  well especially with Catch2 due to the incomplete support for C++17.
#  If neither is chosen, the build relies on the system default, which is
#  on LC with `USE_GCC_LIBCXX` on and `GCC_TOOLCHAIN_VER` set to "4.9.3".
#  If both are on, `USE_GCC_LIBCXX` is turned off. When `USE_GCC_LIBCXX`
#  is on, `GCC_TOOLCHAIN_VER` can be set accordingly (e.g., "8.1.0").

option (USE_CLANG_LIBCXX OFF)
option (USE_GCC_LIBCXX ON)

if (USE_CLANG_LIBCXX)
    set (USE_GCC_LIBCXX OFF)
endif (USE_CLANG_LIBCXX)

if (USE_GCC_LIBCXX)
    set(GCC_TOOLCHAIN_VER "4.9.3" CACHE STRING
        "GCC_TOOLCHAIN_VER chosen by the user at CMake configure time")
    set_property(CACHE GCC_TOOLCHAIN_VER PROPERTY STRINGS
                 "4.9.3" "6.1.0" "7.1.0" "7.3.0" "8.1.0")
endif ()

################################################################
# Initialize build
################################################################

set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
     "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif("${isSystemDir}" STREQUAL "-1")

# Get installation directories -- these get used in various places;
# best to just make them available
include(GNUInstallDirs)
include(ExternalProject)
include(SetupCXX)
include(SetupProtobuf)

################################################################
# Initialize dependencies
################################################################

# Required dependencies

if (NOT WCS_HAS_PROTOBUF)
  if (BUILD_PROTOBUF)
    set(To_Ignore ${BOOST_ROOT}${METIS_ROOT}${SBML_ROOT}${WCS_WITH_METIS})
    message(STATUS "!!!!!!!!!!!!!!!!!!!!!!! Attention !!!!!!!!!!!!!!!!!!!!!!!!!")
    message(STATUS "Type `make` to build protobuf, which is a required package.")
    message(STATUS "Then, rerun cmake in the exactly the same way to configure ")
    message(STATUS "WCS with protobuf just installed.")
    message(STATUS "Alternatively, use -DPROTOBUF_ROOT:PATH=<path-to-protobuf>")
    message(STATUS "to avoid two-staged building.")
    message(STATUS "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    return ()
  else (BUILD_PROTOBUF)
    message(FATAL_ERROR "protobuf is required!")
  endif (BUILD_PROTOBUF)
endif (NOT WCS_HAS_PROTOBUF)

if (DEFINED BOOST_ROOT)
    message(STATUS "BOOST_ROOT: " ${BOOST_ROOT})
    set(Boost_NO_SYSTEM_PATHS ON)
else ()
    if (DEFINED ENV{BOOST_ROOT})
        message(STATUS "ENV BOOST_ROOT: " $ENV{BOOST_ROOT})
        set(Boost_NO_SYSTEM_PATHS ON)
    endif ()
endif ()

find_package(Boost
    #         HINTS ${BOOST_ROOT} $ENV{BOOST_ROOT}
             REQUIRED COMPONENTS
             graph regex filesystem system program_options)

message(STATUS "Boost_INCLUDE_DIRS: " ${Boost_INCLUDE_DIRS})
message(STATUS "Boost_LIBRARY_DIRS: " ${Boost_LIBRARY_DIRS})

# This shouldn't be here, but is ok for now. This will occasionally be
# part of another TPL's libraries (e.g., MKL), but it's no
# guarantee. There's no harm including it multiple times.
if (NOT DL_LIBRARY)
  find_library(DL_LIBRARY dl DOC "The dynamic loader library.")
  if (DL_LIBRARY)
    message(STATUS "Found dl: ${DL_LIBRARY}")
  else ()
    message(FATAL_ERROR
      "dl library not found! This is a required library.\n"
      "Please add the path to libdl to CMAKE_LIBRARY_PATH.")
  endif (DL_LIBRARY)
endif ()

if (WCS_HAS_STD_FILESYSTEM)
  set(LIB_FILESYSTEM "stdc++fs")
else ()
  set(LIB_FILESYSTEM ${Boost_LIBRARIES})
endif ()

include(SetupIHeap)

# Optional dependencies

if (WCS_WITH_SUNDIALS)
  set(WCS_HAS_SUNDIALS FALSE)
  find_package(Sundials MODULE)
  if (SUNDIALS_FOUND)
    set(WCS_HAS_SUNDIALS TRUE)
  endif (SUNDIALS_FOUND)
endif (WCS_WITH_SUNDIALS)

if (WCS_WITH_SBML)
  set(WCS_HAS_SBML FALSE)
  find_package(SBML MODULE)
  if (SBML_FOUND)
    set(WCS_HAS_SBML TRUE)
  endif (SBML_FOUND)
endif (WCS_WITH_SBML)

if (WCS_WITH_EXPRTK)
  include(SetupExprTk)
endif (WCS_WITH_EXPRTK)

if (WCS_WITH_CEREAL)
  include(SetupCereal)
endif (WCS_WITH_CEREAL)

if (WCS_WITH_VTUNE)
  find_package(VTune MODULE)

  if (VTune_FOUND)
    set(WCS_HAS_VTUNE TRUE)
  else ()
    set(WCS_HAS_VTUNE FALSE)
    set(WCS_WITH_VTUNE OFF)
    message(WARNING
      "Requested WCS_WITH_VTUNE=ON, but VTUNE support not detected. "
      "Support NOT enabled. "
      "Try setting VTUNE_DIR to point to the VTune install prefix "
      "and reconfigure.")
  endif (VTune_FOUND)
endif (WCS_WITH_VTUNE)

if (WCS_WITH_HPCTOOLKIT)
  find_package(HPCToolkit MODULE)

  if (HPCToolkit_FOUND)
    set(WCS_HAS_HPCTOOLKIT TRUE)
  else ()
    set(WCS_HAS_HPCTOOLKIT FALSE)
    set(WCS_WITH_HPCTOOLKIT OFF)
    message(WARNING
      "Requested WCS_WITH_HPCTOOLKIT=ON, but HPCTOOLKIT support not detected. "
      "Support NOT enabled. "
      "Try setting HPCTOOLKIT_DIR to point to the HPCToolkit install prefix "
      "and reconfigure.")
  endif (HPCToolkit_FOUND)
endif (WCS_WITH_HPCTOOLKIT)

if (WCS_WITH_ROSS)
  set(WCS_HAS_ROSS FALSE)
  find_package(ROSS MODULE)
  if (ROSS_FOUND)
    find_package(MPI MODULE)
    if (MPI_FOUND)
      set(WCS_HAS_ROSS TRUE)
    endif()
  endif (ROSS_FOUND)
endif (WCS_WITH_ROSS)

if (OpenMP_CXX_FOUND)
  message(STATUS "WCS_OMP_MODE value (${WCS_OMP_MODE})")
  set(WCS_HAS_OPENMP TRUE)
  if(NOT WCS_OMP_MODE)
    set(WCS_OMP_REACTION_UPDATES ON)
  elseif(WCS_OMP_MODE MATCHES 1)
    set(WCS_OMP_REACTION_UPDATES ON)
  elseif(WCS_OMP_MODE MATCHES 2)
    set(WCS_OMP_REACTION_UPDATES ON)
    set(WCS_OMP_REACTION_REACTANTS ON)
  elseif(WCS_OMP_MODE MATCHES 3)
    set(WCS_OMP_REACTION_UPDATES ON)
    set(WCS_OMP_REACTION_REACTANTS ON)
    set(WCS_OMP_REACTION_PRODUCTS ON)
  elseif(WCS_OMP_MODE MATCHES 4)
    set(WCS_OMP_RUN_PARTITION ON)
  elseif(WCS_OMP_MODE MATCHES 5)
    set(WCS_OMP_REACTION_UPDATES ON)
    set(WCS_OMP_RUN_PARTITION ON)
  elseif()
    message(STATUS "Unknown WCS_OMP_MODE (${WCS_OMP_MODE})")
  endif()

  find_package(NUMA MODULE)
  if (NUMA_FOUND)
    set(WCS_HAS_NUMA TRUE)
  endif (NUMA_FOUND)
endif (OpenMP_CXX_FOUND)

if (WCS_WITH_METIS)
  set(WCS_HAS_METIS FALSE)
  find_package(Metis MODULE)
  if (METIS_FOUND)
    set(WCS_HAS_METIS TRUE)
  endif (METIS_FOUND)
endif (WCS_WITH_METIS)

if (WCS_WITH_UNIT_TESTING)
  include(SetupCatch2)
  include(CTest)

  # Now that Catch2 has been found, start adding the unit tests
  # subdirectory addtion has been moved to the CMakeList.txt of each
  # upper directory of unit test directories
  #add_subdirectory(src/reaction_network/unit_test)
  #add_subdirectory(src/utils/unit_test)

  # Add this one last
  #add_subdirectory(unit_test)
endif (WCS_WITH_UNIT_TESTING)

# Handle the documentation
#add_subdirectory(docs)

################################################################
# Build WCS
################################################################

cmake_policy(SET CMP0079 NEW)

# Add WCS source files
add_subdirectory(src)

# Create the WCS library
add_library(wcs ${WCS_SOURCES} ${WCS_HEADERS})
if (WCS_HAS_PROTOBUF)
  add_dependencies(wcs protobuf::libprotobuf)
  add_dependencies(wcs protobuf::protoc)
  if (WCS_HAS_OPENMP)
    add_dependencies(wcs OpenMP::OpenMP_CXX)
    if (WCS_HAS_NUMA)
      add_dependencies(wcs NUMA::numa)
    endif (WCS_HAS_NUMA)
  endif (WCS_HAS_OPENMP)
endif (WCS_HAS_PROTOBUF)

target_include_directories(wcs PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}>)
target_include_directories(wcs SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})

# Use the IMPORTED targets when possible.

if (WCS_HAS_IHEAP)
  add_dependencies(wcs IHEAP-download)
  target_include_directories(wcs PUBLIC
    $<BUILD_INTERFACE:${IHEAP_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}>)
endif (WCS_HAS_IHEAP)

if (WCS_HAS_SUNDIALS)
  target_link_libraries(wcs PUBLIC SUNDIALS::SUNDIALS)
endif (WCS_HAS_SUNDIALS)

if (WCS_HAS_SBML)
  target_link_libraries(wcs PUBLIC SBML::SBML)
#  target_include_directories(wcs SYSTEM PUBLIC ${SBML_INCLUDE_DIRS})
endif (WCS_HAS_SBML)

if (WCS_HAS_EXPRTK)
  add_dependencies(wcs EXPRTK-download)
  target_include_directories(wcs PUBLIC
    $<BUILD_INTERFACE:${EXPRTK_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}>)
endif (WCS_HAS_EXPRTK)

if (WCS_HAS_CEREAL)
  add_dependencies(wcs CEREAL-download)
  target_include_directories(wcs PUBLIC
    $<BUILD_INTERFACE:${CEREAL_HEADER_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}>)
endif (WCS_HAS_CEREAL)

if (WCS_HAS_VTUNE)
  target_link_libraries(wcs PUBLIC ${VTUNE_STATIC_LIB})
  target_include_directories(wcs SYSTEM PUBLIC ${VTUNE_INCLUDE_PATH})
endif ()

if (WCS_HAS_HPCTOOLKIT)
  target_link_libraries(wcs PUBLIC ${HPCTOOLKIT_DYN_LIB})
  target_include_directories(wcs SYSTEM PUBLIC ${HPCTOOLKIT_INCLUDE_PATH})
endif ()

if (WCS_HAS_METIS)
  target_link_libraries(wcs PUBLIC METIS::METIS)
  target_include_directories(wcs SYSTEM PUBLIC ${METIS_INCLUDE_DIR})
endif (WCS_HAS_METIS)

if (WCS_WITH_UNIT_TESTING)
  add_dependencies(wcs CATCH2)
  target_include_directories(wcs PUBLIC
    $<BUILD_INTERFACE:${CATCH2_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/catch2>)
endif (WCS_WITH_UNIT_TESTING)

if (WCS_HAS_PROTOBUF)
  target_link_libraries(wcs PUBLIC wcs_proto protobuf::libprotobuf)
endif (WCS_HAS_PROTOBUF)

if (WCS_HAS_OPENMP)
  target_link_libraries(wcs PUBLIC OpenMP::OpenMP_CXX)
  if (WCS_HAS_NUMA)
    target_link_libraries(wcs PUBLIC NUMA::numa)
  endif (WCS_HAS_NUMA)
endif (WCS_HAS_OPENMP)

################################################################
# Add executables
################################################################

# RPATH setup
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
#https://stackoverflow.com/questions/54464929/cmake-install-rpath-use-link-path-without-any-effect
message(STATUS "CMAKE_INSTALL_RPATH = ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
unset (LIBRARY_PATH)

# add executable reaction
add_executable( reaction-bin src/reaction.cpp )
target_include_directories(reaction-bin PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}>)

target_link_libraries(reaction-bin PRIVATE wcs ${LIB_FILESYSTEM})
set_target_properties(reaction-bin PROPERTIES OUTPUT_NAME reaction)
set_target_properties(reaction-bin PROPERTIES CMAKE_INSTALL_RPATH
                      "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
list(APPEND WCS_EXEC_TARGETS reaction-bin)

# add executable ssa
add_executable( ssa-bin src/ssa.cpp )
target_include_directories(ssa-bin PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(ssa-bin PRIVATE wcs ${LIB_FILESYSTEM})
set_target_properties(ssa-bin PROPERTIES OUTPUT_NAME ssa)
set_target_properties(ssa-bin PROPERTIES CMAKE_INSTALL_RPATH
                      "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
list(APPEND WCS_EXEC_TARGETS ssa-bin)


# add executable ssa_partition
if (WCS_HAS_METIS)
  add_executable( ssa_partition-bin src/ssa_partition.cpp)
  target_include_directories(ssa_partition-bin PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}>)

  target_include_directories(ssa_partition-bin
      PRIVATE $<TARGET_PROPERTY:wcs,INTERFACE_INCLUDE_DIRECTORIES>
      ${Protobuf_DIR}/include ${Boost_INCLUDE_DIRS})

  target_link_libraries(ssa_partition-bin wcs ${LIB_FILESYSTEM}
                        -L${Boost_LIBRARY_DIRS} boost_program_options)

  set_target_properties(ssa_partition-bin PROPERTIES OUTPUT_NAME ssa_partition)
  set_target_properties(ssa_partition-bin PROPERTIES CMAKE_INSTALL_RPATH
                        "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
  list(APPEND WCS_EXEC_TARGETS ssa_partition-bin)
endif (WCS_HAS_METIS)

# add executable graphml2dot
add_executable( graphml2dot-bin src/utils/graphml2dot.cpp )
set_target_properties(graphml2dot-bin PROPERTIES OUTPUT_NAME graphml2dot)
target_link_libraries(graphml2dot-bin PRIVATE wcs ${LIB_FILESYSTEM})

list(APPEND WCS_EXEC_TARGETS graphml2dot-bin)

# add executable partition
add_executable( partition-bin src/partition.cpp )
target_include_directories(partition-bin PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}>)

if (WCS_HAS_METIS)
  target_include_directories(partition-bin SYSTEM PUBLIC ${METIS_INCLUDE_DIR})
endif (WCS_HAS_METIS)

target_link_libraries(partition-bin PRIVATE wcs ${LIB_FILESYSTEM})
set_target_properties(partition-bin PROPERTIES OUTPUT_NAME partition)
set_target_properties(partition-bin PROPERTIES CMAKE_INSTALL_RPATH
                      "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

list(APPEND WCS_EXEC_TARGETS partition-bin)

# add executable loadSBML
if (WCS_HAS_SBML)
  add_executable( loadSBML-bin src/utils/loadSBML.cpp )
  target_compile_definitions(loadSBML-bin PRIVATE WCS_HAS_CONFIG)
  target_include_directories(loadSBML-bin
    PRIVATE $<TARGET_PROPERTY:wcs,INTERFACE_INCLUDE_DIRECTORIES> ${SBML_INCLUDE_DIR})
  target_link_libraries(loadSBML-bin PRIVATE SBML::SBML)
  set_target_properties(loadSBML-bin PROPERTIES OUTPUT_NAME loadSBML)
  list(APPEND WCS_EXEC_TARGETS loadSBML-bin)

  add_executable( generate_cxx_code-bin src/generate_cxx_code.cpp )
  target_compile_definitions(generate_cxx_code-bin PRIVATE WCS_HAS_CONFIG)
  target_include_directories(generate_cxx_code-bin
    PRIVATE $<TARGET_PROPERTY:wcs,INTERFACE_INCLUDE_DIRECTORIES> ${SBML_INCLUDE_DIR})
  target_link_libraries(generate_cxx_code-bin PRIVATE SBML::SBML wcs ${LIB_FILESYSTEM})
  set_target_properties(generate_cxx_code-bin PROPERTIES OUTPUT_NAME generate_cxx_code)
  list(APPEND WCS_EXEC_TARGETS generate_cxx_code-bin)
endif (WCS_HAS_SBML)

# add executable des
if (WCS_HAS_ROSS)
  add_executable( des-bin src/des.cpp )
  target_include_directories(des-bin SYSTEM PUBLIC ${ROSS_INCLUDE_DIRS})
  target_include_directories(des-bin PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}>)
  target_link_libraries(des-bin PRIVATE ROSS::ROSS MPI::MPI_CXX wcs ${LIB_FILESYSTEM})
  set_target_properties(des-bin PROPERTIES OUTPUT_NAME des)
  list(APPEND WCS_EXEC_TARGETS des-bin)
endif (WCS_HAS_ROSS)

# Build tests
add_executable( t_state_rngen-bin src/utils/unit_tests/t_state_rngen.cpp )
target_include_directories(t_state_rngen-bin PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}>)

target_link_libraries(t_state_rngen-bin PRIVATE wcs ${LIB_FILESYSTEM})
set_target_properties(t_state_rngen-bin PROPERTIES OUTPUT_NAME t_state_rngen)
set_target_properties(t_state_rngen-bin PROPERTIES CMAKE_INSTALL_RPATH
                      "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

list(APPEND WCS_UNIT_TEST_TARGETS t_state_rngen-bin)

# Install the binaries
install(
  TARGETS  ${WCS_EXEC_TARGETS}
  EXPORT   WCSTargets
  RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
  ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

install(
  TARGETS  ${WCS_UNIT_TEST_TARGETS}
#  EXPORT   WCSTargets
  RUNTIME  DESTINATION ${CMAKE_INSTALL_PREFIX}/tests/unit
  ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

# Install
install(FILES tests/simple.graphml
  DESTINATION ${CMAKE_INSTALL_DATADIR})

# It is used in the utils/generate_cxx_code.cpp file
set (WCS_SRC_DIR ${CMAKE_SOURCE_DIR})

# Write the configure file
configure_file(
  "${CMAKE_SOURCE_DIR}/cmake/configure_files/wcs_config.hpp.in"
  "${CMAKE_BINARY_DIR}/wcs_config.hpp"
  @ONLY)
target_compile_definitions(wcs PUBLIC WCS_HAS_CONFIG)

if (TARGET WCS_CXX_FLAGS_werror)
  target_link_libraries(wcs PUBLIC WCS_CXX_FLAGS_werror)
endif ()

target_link_libraries(wcs PUBLIC ${DL_LIBRARY})

target_link_libraries(wcs PUBLIC ${Boost_LIBRARIES})

# Clean things up
include(WCSDebugUtilities)
wcs_remove_default_include_paths_from_all_subtargets(wcs)
# This is to fix a bug with certain compilers interacting poorly with
# NVCC. In particular, the NVCC include line would have things like
# "-isystem=/usr/include" which would cause issues resolving STL
# include paths. Since compilers include "/usr/include" in their
# system include paths anyway (things searched by #include <...>), we
# can safely remove these from the explicit link line.

if (WCS_DEBUG_PRINT_SUBTARGETS)
  wcs_print_all_subtargets(wcs)
endif ()

# Add the rest of the things
#add_subdirectory(tests)
#add_subdirectory(scripts)

################################################################
# Install WCS
################################################################

include(CMakePackageConfigHelpers)

# Write the version file. This is independent of build/install tree.
write_basic_package_version_file(
  WCSConfigVersion.cmake
  VERSION "${WCS_VERSION}"
  COMPATIBILITY SameMajorVersion)

# This is for the build tree
set(INCLUDE_INSTALL_DIRS
  "${CMAKE_SOURCE_DIR}/src"
  "${CMAKE_SOURCE_DIR}/src/params"
  "${CMAKE_SOURCE_DIR}/src/partition"
  "${CMAKE_SOURCE_DIR}/src/proto"
  "${CMAKE_SOURCE_DIR}/src/reaction_network"
  "${CMAKE_SOURCE_DIR}/src/sim_methods"
  "${CMAKE_SOURCE_DIR}/src/utils")
set(LIB_INSTALL_DIR "${CMAKE_BINARY_DIR}")
set(EXTRA_CMAKE_MODULE_DIR "${CMAKE_SOURCE_DIR}/cmake/modules")
configure_package_config_file(cmake/configure_files/WCSConfig.cmake.in
  "${CMAKE_BINARY_DIR}/WCSConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_BINARY_DIR}"
  PATH_VARS INCLUDE_INSTALL_DIRS LIB_INSTALL_DIR)

# Write the configure file for the install tree
set(INCLUDE_INSTALL_DIRS include)
set(LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})
set(CMAKE_INSTALL_DIR ${LIB_INSTALL_DIR}/cmake/wcs)
set(EXTRA_CMAKE_MODULE_DIR)
configure_package_config_file(cmake/configure_files/WCSConfig.cmake.in
  "${CMAKE_BINARY_DIR}/WCSConfig.cmake.install"
  INSTALL_DESTINATION "${CMAKE_INSTALL_DIR}"
  PATH_VARS INCLUDE_INSTALL_DIRS LIB_INSTALL_DIR)


# Install library
install(
  TARGETS wcs
  EXPORT WCSTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )

# Build tree export
#export(EXPORT WCSTargets NAMESPACE WCS:: FILE WCSTargets.cmake)

# Install export
install(EXPORT WCSTargets
  NAMESPACE WCS::
  FILE WCSTargets.cmake
  DESTINATION "${CMAKE_INSTALL_DIR}")

# Install the cmake stuff
install(FILES
  "${PROJECT_BINARY_DIR}/WCSConfig.cmake.install"
  DESTINATION "${CMAKE_INSTALL_DIR}"
  RENAME "WCSConfig.cmake")
install(FILES "${PROJECT_BINARY_DIR}/WCSConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_DIR})
install(DIRECTORY cmake/modules
  DESTINATION "${CMAKE_INSTALL_DIR}"
  FILES_MATCHING PATTERN "*.cmake")

# Install header files
install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/src/params"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/src/partition"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/src/proto"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/src/reaction_network"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/src/sim_methods"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/src/utils"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
install(
  FILES "${PROJECT_BINARY_DIR}/wcs_config.hpp"
        "${PROJECT_SOURCE_DIR}/src/wcs_types.hpp"
        "${PROJECT_SOURCE_DIR}/src/bgl.hpp"
        "${PROJECT_SOURCE_DIR}/src/des.hpp"
        "${PROJECT_SOURCE_DIR}/src/wcs-ross-bf.hpp"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

if (WCS_HAS_EXPRTK)
  install(
    FILES ${EXPRTK_DIR}/${EXPRTK_HEADER}
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
endif ()

# TODO: either need to install the whole directory or not at all
if (WCS_HAS_CEREAL)
  install(
    DIRECTORY "${CEREAL_HEADER_DIR}/cereal"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
endif ()

# Install contributor list, license, readme
install(
  FILES #"${PROJECT_SOURCE_DIR}/CONTRIBUTORS"
        #"${PROJECT_SOURCE_DIR}/LICENSE"
        "${PROJECT_SOURCE_DIR}/README.md"
  DESTINATION "${CMAKE_INSTALL_DOCDIR}")

################################################################
# Configuration summary
################################################################

# This creates a formatted string that contains a list of variables,
# one per line, with their values interpreted as TRUE or FALSE. The
# purpose is to provide uniform output, rather than an odd mixture of
# "1", "0", "ON", "OFF", "TRUE" and "FALSE".
macro(append_str_tf STRING_VAR)
  wcs_get_max_str_length(_max_length ${ARGN})
  math(EXPR _max_length "${_max_length} + 2")

  foreach(var ${ARGN})
    string(LENGTH "${var}" _var_length)
    math(EXPR _num_spaces "${_max_length} - ${_var_length}")
    wcs_get_space_string(_spaces ${_num_spaces})
    if (${var})
      set(${var} "TRUE")
      string(APPEND ${STRING_VAR} "  ${var}:" "${_spaces}" "TRUE\n")
    else ()
      set(${var} "FALSE")
      string(APPEND ${STRING_VAR} "  ${var}:" "${_spaces}" "FALSE\n")
    endif ()
  endforeach()
endmacro ()

# NOTE: message() outputs to stderr by default. We now use a string to
# maintain this information and then have cmake echo it to stdout. The
# only side effects are that if you use the CMake GUI, you won't see
# this output anymore (they only report stderr) and that if you add
# something to the list, you must remember your newline!
set(_str "\n== WCS Configuration Summary ==\n\n")
string(APPEND _str "  PROJECT_SOURCE_DIR:   ${PROJECT_SOURCE_DIR}\n"
  "  PROJECT_BINARY_DIR:   ${PROJECT_BINARY_DIR}\n\n"
  "  CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}\n"
  "  CMAKE_BUILD_TYPE:     ${CMAKE_BUILD_TYPE}\n\n")
if (CMAKE_BUILD_TYPE MATCHES None)
  string(APPEND _str
    "  CXX FLAGS:            ${CMAKE_CXX_FLAGS}\n")
elseif (CMAKE_BUILD_TYPE MATCHES Release)
  string(APPEND _str
    "  CXX FLAGS:            ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}\n")
elseif (CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
  string(APPEND _str
    "  CXX FLAGS:            ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}\n")
elseif (CMAKE_BUILD_TYPE MATCHES Debug)
  string(APPEND _str
    "  CXX FLAGS:            ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}\n")
endif ()
string(APPEND _str "\n")

#Print the true/false guys
append_str_tf(_str
  WCS_GNU_LINUX
  WCS_64BIT_CNT
  WCS_CACHE_DEPENDENT
  WCS_HAS_SUNDIALS
  WCS_HAS_SBML
  WCS_HAS_EXPRTK
  WCS_HAS_CEREAL
  WCS_HAS_CATCH2
  WCS_HAS_VTUNE
  WCS_HAS_HPCTOOLKIT
  WCS_HAS_DOXYGEN
  WCS_HAS_ROSS
  WCS_HAS_OPENMP
  WCS_HAS_NUMA
  WCS_HAS_METIS
  WCS_HAS_PROTOBUF)
string(APPEND _str
  "\n== End WCS Configuration Summary ==\n")

# Output to stdout
execute_process(COMMAND ${CMAKE_COMMAND} -E echo "${_str}")
set(_str)

#
# Write a basic modulefile
#
set(WCS_MODULEFILE_NAME "wcs-${WCS_VERSION}.lua"
  CACHE STRING
  "The name of the WCS modulefile to install. Must end in .lua.")

if (NOT (WCS_MODULEFILE_NAME MATCHES ".+\.lua"))
  message(WARNING
    "WCS_MODULEFILE_NAME must have extension \".lua\". Appending.")
  set(WCS_MODULEFILE_NAME "${WCS_MODULEFILE_NAME}.lua"
    CACHE STRING "" FORCE)
endif ()

configure_file(
  "${CMAKE_SOURCE_DIR}/cmake/configure_files/wcs_module.lua.in"
  "${CMAKE_BINARY_DIR}/wcs_module.lua.install"
  @ONLY)
install(FILES "${CMAKE_BINARY_DIR}/wcs_module.lua.install"
  RENAME "${WCS_MODULEFILE_NAME}"
  DESTINATION "${CMAKE_INSTALL_SYSCONFDIR}/modulefiles")
