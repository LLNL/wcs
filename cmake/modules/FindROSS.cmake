# Defines the following variables:
#   - ROSS_FOUND
#   - ROSS_LIBRARY
#   - ROSS_INCLUDE_DIR
#
# Also creates an imported target ROSS

# Find the header
find_path(ROSS_HEADER_DIR ross.h
  HINTS ${ROSS_ROOT} $ENV{ROSS_ROOT}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with ROSS header.")
find_path(ROSS_HEADER_DIR ross.h)

unset(ROSS_INCLUDE_DIR CACHE)
get_filename_component(ROSS_INCLUDE_DIR ${ROSS_HEADER_DIR}/include DIRECTORY CACHE)

message(STATUS "ROSS_INCLUDE_DIR: ${ROSS_INCLUDE_DIR}")
# Find the library
find_library(ROSS_LIBRARY ROSS
  HINTS ${ROSS_ROOT} $ENV{ROSS_ROOT}
  PATH_SUFFIXES lib64 lib
  NO_DEFAULT_PATH
  DOC "The ROSS library.")
find_library(ROSS_LIBRARY ROSS)

# Standard handling of the package arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ROSS
  DEFAULT_MSG
  ROSS_LIBRARY ROSS_INCLUDE_DIR)

# Setup the imported target
if (NOT TARGET ROSS::ROSS)
  add_library(ROSS::ROSS INTERFACE IMPORTED)
endif (NOT TARGET ROSS::ROSS)

# Set the include directories for the target
set_property(TARGET ROSS::ROSS APPEND
  PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${ROSS_INCLUDE_DIR})

# Set the link libraries for the target
set_property(TARGET ROSS::ROSS APPEND
  PROPERTY INTERFACE_LINK_LIBRARIES ${ROSS_LIBRARY})

#
# Cleanup
#

# Set the include directories
mark_as_advanced(FORCE ROSS_INCLUDE_DIR)

# Set the libraries
set(ROSS_LIBRARIES ROSS::ROSS)
mark_as_advanced(FORCE ROSS_LIBRARY)
