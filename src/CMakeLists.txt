# Add the source files for this directory
set_full_path(THIS_DIR_HEADERS
  bgl.hpp
  wcs_types.hpp
  des.hpp
  ssa-cfg.hpp
  wcs-ross-bf.hpp
  wcs_params.hpp
  )

set_full_path(THIS_DIR_SOURCES
  wcs_params.cpp
  )

# Add the subdirectories
add_subdirectory(partition)
add_subdirectory(proto)
add_subdirectory(reaction_network)
add_subdirectory(sim_methods)
add_subdirectory(utils)

# Propagate the files up the tree
set(WCS_HEADERS "${WCS_HEADERS}" "${THIS_DIR_HEADERS}" PARENT_SCOPE)
set(WCS_SOURCES "${WCS_SOURCES}" "${THIS_DIR_SOURCES}" PARENT_SCOPE)
