# Add the header and source files for this directory
set_full_path(THIS_DIR_HEADERS
  detect_methods.hpp
  exception.hpp
  graph_factory.hpp
  write_graphviz.hpp
  )

set_full_path(THIS_DIR_SOURCES
  exception.cpp
  graph_factory.cpp
  graphml2dot.cpp
  )

# Add the subdirectories

# Propagate the files up the tree
set(WCS_HEADERS "${HEADERS}" "${THIS_DIR_HEADERS}" PARENT_SCOPE)
set(WCS_SOURCES "${SOURCES}" "${THIS_DIR_SOURCES}" PARENT_SCOPE)


add_executable( graphml2dot-bin graphml2dot.cpp )
target_link_libraries(graphml2dot-bin graphml2dot )
set_target_properties(graphml2dot-bin PROPERTIES OUTPUT_NAME graphml2dot)

# Install the binaries
install(
  TARGETS  graphml2dot-bin
  EXPORT   WCSTargets
  RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
  ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

# Install the relevant prototext
#install(DIRECTORY unit_tests/examples 
#  DESTINATION ${CMAKE_INSTALL_DATADIR}/examples)
