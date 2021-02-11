# Add the header and source files for this directory
set_full_path(THIS_DIR_HEADERS
  ssa_params.hpp
  wcs_params.hpp
  )

set_full_path(THIS_DIR_SOURCES
  ssa_params.cpp
  wcs_params.cpp
  )

# Add the subdirectories

# Propagate the files up the tree
set(WCS_HEADERS "${WCS_HEADERS}" "${THIS_DIR_HEADERS}" PARENT_SCOPE)
set(WCS_SOURCES "${WCS_SOURCES}" "${THIS_DIR_SOURCES}" PARENT_SCOPE)
