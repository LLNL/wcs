# Add the header and source files for this directory
set_full_path(THIS_DIR_HEADERS
  detect_methods.hpp
  exception.hpp
  graph_factory.hpp
  rngen.hpp
  rngen_impl.hpp
  trace.hpp
  trace_ssa.hpp
  samples.hpp
  write_graphviz.hpp
  streamvec.hpp
  streamvec_impl.hpp
  streambuff.hpp
  streambuff_impl.hpp
  state_io.hpp
  state_io_impl.hpp
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


add_executable( graphml2dot-bin graphml2dot.cpp )
target_link_libraries(graphml2dot-bin PRIVATE wcs)
set_target_properties(graphml2dot-bin PROPERTIES OUTPUT_NAME graphml2dot)

list(APPEND WCS_UTILS_TARGETS graphml2dot-bin)

if (WCS_HAS_SBML)
  add_executable( readSBML-bin readSBML.cpp )
  target_compile_definitions(readSBML-bin PRIVATE WCS_HAS_CONFIG)
  target_include_directories(readSBML-bin
    PRIVATE $<TARGET_PROPERTY:wcs,INTERFACE_INCLUDE_DIRECTORIES>)
  target_link_libraries(readSBML-bin PRIVATE SBML::SBML)
  set_target_properties(readSBML-bin PROPERTIES OUTPUT_NAME readSBML)
  list(APPEND WCS_UTILS_TARGETS readSBML-bin)
endif (WCS_HAS_SBML)

# Install the binaries
install(
  TARGETS  ${WCS_UTILS_TARGETS}
  EXPORT   WCSTargets
  RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
  ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

# Install the relevant prototext
#install(DIRECTORY unit_tests/examples
#  DESTINATION ${CMAKE_INSTALL_DATADIR}/examples)
