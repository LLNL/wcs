# Defines the following variables:
#   - SUNDIALS_FOUND
#   - SUNDIALS_LIBRARY
#   - SUNDIALS_INCLUDE_DIR
#
# Also creates an imported target SUNDIALS

# Find the header
find_path(SUNDIALS_HEADER_DIR cvode.h
  HINTS ${SUNDIALS_ROOT} $ENV{SUNDIALS_ROOT}
  PATH_SUFFIXES include/cvode
  NO_DEFAULT_PATH
  DOC "Directory with Sundials CVODE header.")
find_path(SUNDIALS_HEADER_DIR cvode.h)

unset(SUNDIALS_INCLUDE_DIR CACHE)
get_filename_component(SUNDIALS_INCLUDE_DIR ${SUNDIALS_HEADER_DIR} DIRECTORY CACHE)

message(STATUS "SUNDIALS_INCLUDE_DIR: ${SUNDIALS_INCLUDE_DIR}")
# Find the library
find_library(SUNDIALS_LIBRARY sundials_cvode
  HINTS ${SUNDIALS_ROOT} $ENV{SUNDIALS_ROOT}
  PATH_SUFFIXES lib64 lib
  NO_DEFAULT_PATH
  DOC "The Sundials CVODE library.")
find_library(SUNDIALS_LIBRARY sundials_cvode)

# Standard handling of the package arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SUNDIALS
  DEFAULT_MSG
  SUNDIALS_LIBRARY SUNDIALS_INCLUDE_DIR)

# Setup the imported target
if (NOT TARGET SUNDIALS::SUNDIALS)
  add_library(SUNDIALS::SUNDIALS INTERFACE IMPORTED)
endif (NOT TARGET SUNDIALS::SUNDIALS)

# Set the include directories for the target
set_property(TARGET SUNDIALS::SUNDIALS APPEND
  PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${SUNDIALS_INCLUDE_DIR})

# Set the link libraries for the target
set_property(TARGET SUNDIALS::SUNDIALS APPEND
  PROPERTY INTERFACE_LINK_LIBRARIES ${SUNDIALS_LIBRARY})

#
# Cleanup
#

# Set the include directories
mark_as_advanced(FORCE SUNDIALS_INCLUDE_DIR)

# Set the libraries
set(SUNDIALS_LIBRARIES SUNDIALS::SUNDIALS)
mark_as_advanced(FORCE SUNDIALS_LIBRARY)
