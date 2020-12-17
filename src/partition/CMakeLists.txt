# Add the header and source files for this directory
set_full_path(THIS_DIR_HEADERS
  metis_params.hpp
  metis_partition.hpp
  partition.hpp
  partition_info.hpp
  )

set_full_path(THIS_DIR_SOURCES
  metis_params.cpp
  metis_partition.cpp
  partition.cpp
  partition_info.cpp
  )

# Add the subdirectories

# Propagate the files up the tree
set(WCS_HEADERS "${WCS_HEADERS}" "${THIS_DIR_HEADERS}" PARENT_SCOPE)
set(WCS_SOURCES "${WCS_SOURCES}" "${THIS_DIR_SOURCES}" PARENT_SCOPE)
