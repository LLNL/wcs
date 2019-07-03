# Add the source files for this directory
set_full_path(THIS_DIR_HEADERS
  bgl.hpp
  wcs_types.hpp
  )

set_full_path(THIS_DIR_SOURCES
  )

# Add the subdirectories
add_subdirectory(reaction_network)
add_subdirectory(sim_methods)
add_subdirectory(utils)

# Propagate the files up the tree
set(WCS_HEADERS "${WCS_HEADERS}" "${THIS_DIR_HEADERS}" PARENT_SCOPE)
set(WCS_SOURCES "${WCS_SOURCES}" "${THIS_DIR_SOURCES}" PARENT_SCOPE)

# Add executable
add_executable( reaction-bin reaction.cpp )
target_include_directories(reaction-bin PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}>)

target_link_libraries(reaction-bin PRIVATE wcs)
set_target_properties(reaction-bin PROPERTIES OUTPUT_NAME reaction)

add_executable( ssa-bin ssa.cpp )
target_include_directories(ssa-bin PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}>)

target_link_libraries(ssa-bin PRIVATE wcs)
set_target_properties(ssa-bin PROPERTIES OUTPUT_NAME ssa)

# Install the binaries
install(
  TARGETS  reaction-bin ssa-bin
  EXPORT   WCSTargets
  RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
  ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

# Install
install(DIRECTORY utils/unit_tests/examples
  DESTINATION ${CMAKE_INSTALL_DATADIR}/examples)
