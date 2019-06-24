# Defines the following variables:
#   - SBML_FOUND
#   - SBML_LIBRARIES
#   - SBML_INCLUDE_DIRS
#
# Also creates an imported target SBML

# Find the header
find_path(SBML_INCLUDE_DIRS sbml.h
  HINTS ${SBML_DIR} $ENV{SBML_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with SBML header.")
find_path(SBML_INCLUDE_DIRS sbml.h)

# Find the library
find_library(SBML_LIBRARY sbml
  HINTS ${SBML_DIR} $ENV{SBML_DIR}
  PATH_SUFFIXES lib64 lib
  NO_DEFAULT_PATH
  DOC "The SBML library.")
find_library(SBML_LIBRARY sbml)

# Standard handling of the package arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SBML
  DEFAULT_MSG
  SBML_LIBRARY SBML_INCLUDE_DIRS)

# Setup the imported target
if (NOT TARGET SBML::SBML)
  add_library(SBML::SBML INTERFACE IMPORTED)
endif (NOT TARGET SBML::SBML)

# Set the include directories for the target
set_property(TARGET SBML::SBML APPEND
  PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${SBML_INCLUDE_DIRS})

# Set the link libraries for the target
set_property(TARGET SBML::SBML APPEND
  PROPERTY INTERFACE_LINK_LIBRARIES ${SBML_LIBRARY})

#
# Cleanup
#

# Set the include directories
mark_as_advanced(FORCE SBML_INCLUDE_DIRS)

# Set the libraries
set(SBML_LIBRARIES SBML::SBML)
mark_as_advanced(FORCE SBML_LIBRARY)
