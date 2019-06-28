# Add the header and source files for this directory
set_full_path(THIS_DIR_HEADERS
  ssa_nrm.hpp
  )

set_full_path(THIS_DIR_SOURCES
  ssa_nrm.cpp
  )

# Propagate the files up the tree
set(WCS_HEADERS "${HEADERS}" "${THIS_DIR_HEADERS}" PARENT_SCOPE)
set(WCS_SOURCES "${SOURCES}" "${THIS_DIR_SOURCES}" PARENT_SCOPE)
