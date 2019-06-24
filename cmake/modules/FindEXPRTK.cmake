# Defines the following variables:
#   - EXPRTK_FOUND
#   - EXPRTK_LIBRARIES
#   - EXPRTK_INCLUDE_DIRS
#
# Also creates an imported target EXPRTK

# Find the header
find_path(EXPRTK_INCLUDE_DIRS exprtk.h
  HINTS ${EXPRTK_DIR} $ENV{EXPRTK_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with EXPRTK header.")
find_path(EXPRTK_INCLUDE_DIRS exprtk.h)

# Find the library
find_library(EXPRTK_LIBRARY exprtk
  HINTS ${EXPRTK_DIR} $ENV{EXPRTK_DIR}
  PATH_SUFFIXES lib64 lib
  NO_DEFAULT_PATH
  DOC "The EXPRTK library.")
find_library(EXPRTK_LIBRARY exprtk)

# Standard handling of the package arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EXPRTK
  DEFAULT_MSG
  EXPRTK_LIBRARY EXPRTK_INCLUDE_DIRS)

# Setup the imported target
if (NOT TARGET EXPRTK::EXPRTK)
  add_library(EXPRTK::EXPRTK INTERFACE IMPORTED)
endif (NOT TARGET EXPRTK::EXPRTK)

# Set the include directories for the target
set_property(TARGET EXPRTK::EXPRTK APPEND
  PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${EXPRTK_INCLUDE_DIRS})

# Set the link libraries for the target
set_property(TARGET EXPRTK::EXPRTK APPEND
  PROPERTY INTERFACE_LINK_LIBRARIES ${EXPRTK_LIBRARY})

#
# Cleanup
#

# Set the include directories
mark_as_advanced(FORCE EXPRTK_INCLUDE_DIRS)

# Set the libraries
set(EXPRTK_LIBRARIES EXPRTK::EXPRTK)
mark_as_advanced(FORCE EXPRTK_LIBRARY)
