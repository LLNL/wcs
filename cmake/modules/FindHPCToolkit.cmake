# Exports the following variables
#   HPCToolkit_FOUND
#   HPCTOOLKIT_INCLUDE_DIR
#   HPCTOOLKIT_DYN_LIB
#

find_path(HPCToolkit_INCLUDE_PATH hpctoolkit.h
  HINTS ${HPCToolkit_DIR} $ENV{HPCToolkit_DIR}
        ${HPCTOOLKIT_DIR} $ENV{HPCTOOLKIT_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "The location of HPCToolkit headers.")
find_path(HPCToolkit_INCLUDE_PATH libittnotify.h)

find_library(HPCTOOLKIT_DYN_LIB libhpctoolkit.so
  HINTS ${HPCToolkit_DIR} $ENV{HPCToolkit_DIR}
        ${HPCTOOLKIT_DIR} $ENV{HPCTOOLKIT_DIR}
  PATH_SUFFIXES lib/hpctoolkit
  NO_DEFAULT_PATH
  DOC "The location of HPCToolkit library.")
find_library(HPCTOOLKIT_DYN_LIB libhpctoolkit.so)

# Standard handling of the package arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HPCToolkit
  REQUIRED_VARS HPCToolkit_INCLUDE_PATH HPCTOOLKIT_DYN_LIB)

if (HPCToolkit_INCLUDE_PATH AND HPCTOOLKIT_DYN_LIB)
  if (NOT TARGET hpctoolkit::hpctoolkit)
    add_library(hpctoolkit::hpctoolkit SHARED IMPORTED)
  endif ()

  set_target_properties(hpctoolkit::hpctoolkit PROPERTIES
    IMPORTED_LOCATION "${HPCTOOLKIT_DYN_LIB}"
    INTERFACE_INCLUDE_DIRECTORIES "${HPCToolkit_INCLUDE_PATH}")

  set(HPCToolkit_LIBRARIES hpctoolkit::hpctoolkit)
endif ()
