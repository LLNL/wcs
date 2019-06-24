# Add the header and source files for this directory
set_full_path(THIS_DIR_HEADERS
  vertex_property_base.hpp
  species.hpp
  reaction_base.hpp
  reaction.hpp
  reaction_impl.hpp
  vertex.hpp
  vertex_flat.hpp
  edge.hpp
  network.hpp
  )

set_full_path(THIS_DIR_SOURCES
  vertex_property_base.cpp
  species.cpp
  reaction_base.cpp
  vertex.cpp
  vertex_flat.cpp
  edge.cpp
  network.cpp
  )

# Propagate the files up the tree
set(WCS_HEADERS "${HEADERS}" "${THIS_DIR_HEADERS}" PARENT_SCOPE)
set(WCS_SOURCES "${SOURCES}" "${THIS_DIR_SOURCES}" PARENT_SCOPE)
