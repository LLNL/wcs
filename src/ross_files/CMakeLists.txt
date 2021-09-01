# blt_add_executable(
#     NAME ross_files
#     SOURCES actor.cc
#     DEPENDS_ON ROSS mpi
# )
# find_package(Python COMPONENTS Interpreter)

set(actor_srcs ${CMAKE_CURRENT_SOURCE_DIR}/Funnel.hh)

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/actor.gen.cc
                   COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/generateInfrastructure.py ${actor_srcs} > ${CMAKE_CURRENT_BINARY_DIR}/actor.gen.cc
                   MAIN_DEPENDENCY generateInfrastructure.py
                   DEPENDS generateInfrastructure.py ${actor_srcs}
                   VERBATIM)

set_source_files_properties(
  ${CMAKE_CURRENT_BINARY_DIR}/actor.gen.cc
  PROPERTIES GENERATED TRUE
)

if (WCS_HAS_ROSS)
  add_executable(hybrid-bin
              hybrid.cpp
              hybrid.hpp
              actor.hh
              Checkpoint.hh
              simtime.hh
              actor.cc
              ${actor_srcs}
              ${CMAKE_CURRENT_BINARY_DIR}/actor.gen.cc
              )
  target_link_libraries(hybrid-bin PUBLIC ROSS::ROSS MPI::MPI_CXX wcs ${LIB_FILESYSTEM})
  target_include_directories(hybrid-bin SYSTEM PUBLIC ${ROSS_INCLUDE_DIRS})
  target_include_directories(hybrid-bin PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}>)

  set_target_properties(hybrid-bin PROPERTIES OUTPUT_NAME hybrid)
endif()

set_full_path(THIS_DIR_HEADERS
  )

set_full_path(THIS_DIR_SOURCES 
  )

# # Add the subdirectories


# set(THIS_DIR_SOURCES "${THIS_DIR_SOURCES}" "${CMAKE_CURRENT_BINARY_DIR}/actor.gen.cc")  

# # Propagate the files up the tree
set(WCS_HEADERS "${WCS_HEADERS}" "${THIS_DIR_HEADERS}" PARENT_SCOPE)
set(WCS_SOURCES "${WCS_SOURCES}" "${THIS_DIR_SOURCES}" PARENT_SCOPE)


