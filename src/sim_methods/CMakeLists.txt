# Add the header and source files for this directory
set_full_path(THIS_DIR_HEADERS
  sim_method.hpp
  ssa_nrm.hpp
  )

set_full_path(THIS_DIR_SOURCES
  sim_method.cpp
  ssa_nrm.cpp
  )

# Propagate the files up the tree
set(WCS_HEADERS "${WCS_HEADERS}" "${THIS_DIR_HEADERS}" PARENT_SCOPE)
set(WCS_SOURCES "${WCS_SOURCES}" "${THIS_DIR_SOURCES}" PARENT_SCOPE)
