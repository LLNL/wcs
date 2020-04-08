# Add the header and source files for this directory
set_full_path(THIS_DIR_HEADERS
  detect_methods.hpp
  exception.hpp
  graph_factory.hpp
  rngen.hpp
  rngen_impl.hpp
  timer.hpp
  trace.hpp
  trace_ssa.hpp
  samples.hpp
  seed.hpp
  write_graphviz.hpp
  streamvec.hpp
  streamvec_impl.hpp
  streambuff.hpp
  streambuff_impl.hpp
  state_io.hpp
  state_io_impl.hpp
  state_io_cereal.hpp
  )

set_full_path(THIS_DIR_SOURCES
  exception.cpp
  graph_factory.cpp
  trace.cpp
  trace_ssa.cpp
  samples.cpp
  )

# Add the subdirectories

# Propagate the files up the tree
set(WCS_HEADERS "${WCS_HEADERS}" "${THIS_DIR_HEADERS}" PARENT_SCOPE)
set(WCS_SOURCES "${WCS_SOURCES}" "${THIS_DIR_SOURCES}" PARENT_SCOPE)
