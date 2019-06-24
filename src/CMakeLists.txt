# Add the source files for this directory
set_full_path(THIS_DIR_HEADERS
  bgl.hpp
  wcs_types.hpp
  )

set_full_path(THIS_DIR_SOURCES
  reaction.cpp
  )

# Add the subdirectories
add_subdirectory(reaction_network)
add_subdirectory(utils)

# Propagate the files up the tree
set(WCS_HEADERS "${HEADERS}" "${THIS_DIR_HEADERS}" PARENT_SCOPE)
set(WCS_SOURCES "${SOURCES}" "${THIS_DIR_SOURCES}" PARENT_SCOPE)

add_executable( reaction-bin reaction.cpp )
target_link_libraries(reaction-bin reaction )
set_target_properties(reaction-bin PROPERTIES OUTPUT_NAME reaction)

# Install the binaries
install(
  TARGETS  reaction-bin
  EXPORT   WCSTargets
  RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
  ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

# Install the relevant prototext
install(DIRECTORY utils/unit_tests/examples 
  DESTINATION ${CMAKE_INSTALL_DATADIR}/examples)
