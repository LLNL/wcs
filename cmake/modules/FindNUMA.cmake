# Exports the following variables
#   NUMA_FOUND
#   NUMA_INCLUDE_PATH
#   NUMA_LIBRARY
#

find_path(NUMA_INCLUDE_PATH numa.h
  HINTS ${NUMA_DIR} $ENV{NUMA_DIR}
        ${NUMA_ROOT} $ENV{NUMA_ROOT}
        "/usr" "/"
  PATH_SUFFIXES include
  DOC "The location of NUMA headers.")
find_path(NUMA_INCLUDE_PATH numa.h)
message(STATUS "NUMA_INCLUDE_PATH: ${NUMA_INCLUDE_PATH}")

find_library(NUMA_LIBRARY numa
  HINTS ${NUMA_DIR} $ENV{NUMA_DIR}
        ${NUMA_ROOT} $ENV{NUMA_ROOT}
        "/usr" "/"
  PATH_SUFFIXES lib64 lib
  DOC "The location of NUMA Static library.")
find_library(NUMA_LIBRARY numa)
message(STATUS "NUMA_LIBRARY: ${NUMA_LIBRARY}")

# Standard handling of the package arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NUMA
  REQUIRED_VARS NUMA_INCLUDE_PATH NUMA_LIBRARY)

if (NUMA_INCLUDE_PATH AND NUMA_LIBRARY)
  if (NOT TARGET NUMA::numa)
    add_library(NUMA::numa SHARED IMPORTED)
  endif ()

  set_target_properties(NUMA::numa PROPERTIES
    IMPORTED_LOCATION "${NUMA_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${NUMA_INCLUDE_PATH}")

  set(NUMA_LIBRARIES NUMA::numa)
endif ()
