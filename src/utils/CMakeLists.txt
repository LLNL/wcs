# Add the header and source files for this directory
set_full_path(THIS_DIR_HEADERS
  detect_methods.hpp
  exception.hpp
  file.hpp
  generate_cxx_code.hpp
  graph_factory.hpp
  input_filetype.hpp
  omp_diagnostics.hpp
  print_vertices.hpp
  rngen.hpp
  rngen_impl.hpp
  samples_ssa.hpp
  sbml_utils.hpp
  sbml_functions.hpp
  seed.hpp
  state_io.hpp
  state_io_impl.hpp
  state_io_cereal.hpp
  streambuff.hpp
  streambuff_impl.hpp
  streamvec.hpp
  streamvec_impl.hpp
  timer.hpp
  trajectory.hpp
  trace_ssa.hpp
  trace_generic.hpp
  traits.hpp
  write_graphviz.hpp
  write_graphviz_impl.hpp
  )

set_full_path(THIS_DIR_SOURCES
  exception.cpp
  file.cpp
  generate_cxx_code.cpp
  graph_factory.cpp
  input_filetype.cpp
  omp_diagnostics.cpp
  print_vertices.cpp
  samples_ssa.cpp
  sbml_utils.cpp
  sbml_functions.cpp
  trajectory.cpp
  trace_ssa.cpp
  trace_generic.cpp
  )

# Add the subdirectories
#if (WCS_WITH_UNIT_TESTING)
  add_subdirectory(unit_tests)
#endif ()

# Propagate the files up the tree
set(WCS_HEADERS "${WCS_HEADERS}" "${THIS_DIR_HEADERS}" PARENT_SCOPE)
set(WCS_SOURCES "${WCS_SOURCES}" "${THIS_DIR_SOURCES}" PARENT_SCOPE)
