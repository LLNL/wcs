# Defines the following variables:
#   - METIS_FOUND
#   - METIS_LIBRARY
#   - METIS_INCLUDE_DIR
#
# Also creates an imported target METIS

# Find the header
find_path(METIS_HEADER_DIR metis.h
  HINTS ${METIS_ROOT}/include $ENV{METIS_ROOT}/include
  NO_DEFAULT_PATH
  DOC "Directory with METIS header.")
find_path(METIS_HEADER_DIR metis.h)

unset(METIS_INCLUDE_DIR CACHE)
get_filename_component(METIS_INCLUDE_DIR ${METIS_HEADER_DIR}/include DIRECTORY CACHE)

message(STATUS "METIS_INCLUDE_DIR: ${METIS_INCLUDE_DIR}")
# Find the library
find_library(METIS_LIBRARY metis
  HINTS ${METIS_ROOT} $ENV{METIS_ROOT}
  PATH_SUFFIXES lib64 lib
  NO_DEFAULT_PATH
  DOC "The METIS library.")
find_library(METIS_LIBRARY metis)

# Standard handling of the package arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Metis
  DEFAULT_MSG
  METIS_LIBRARY METIS_INCLUDE_DIR)

# Setup the imported target
if (NOT TARGET METIS::METIS)
  add_library(METIS::METIS INTERFACE IMPORTED)
endif (NOT TARGET METIS::METIS)

# Set the include directories for the target
set_property(TARGET METIS::METIS APPEND
  PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${METIS_INCLUDE_DIR})

# Set the link libraries for the target
set_property(TARGET METIS::METIS APPEND
  PROPERTY INTERFACE_LINK_LIBRARIES ${METIS_LIBRARY})

#
# Cleanup
#

# Set the include directories
mark_as_advanced(FORCE METIS_INCLUDE_DIR)

# Set the libraries
set(METIS_LIBRARIES METIS::METIS)
mark_as_advanced(FORCE METIS_LIBRARY)
