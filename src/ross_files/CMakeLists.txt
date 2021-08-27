# blt_add_executable(
#     NAME ross_files
#     SOURCES actor.cc
#     DEPENDS_ON ROSS mpi
# )
# find_package(Python COMPONENTS Interpreter)

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/actor.gen.cc
                   COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/generateInfrastructure.py ${actor_srcs} > ${CMAKE_CURRENT_BINARY_DIR}/actor.gen.cc
                   MAIN_DEPENDENCY generateInfrastructure.py
                   DEPENDS generateInfrastructure.py ${actor_srcs}
                   VERBATIM)
SET_PROPERTY(SOURCE ${CMAKE_CURRENT_BINARY_DIR}/actor.gen.cc
             PROPERTY INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR})


# SET(my_srcs
#     ${actor_srcs}
#     actor.hh
#     actor.cc
#     simtime.hh
#     Checkpoint.hh
# )

# # copy .hh and .cc fiels in the folder with actor.gen.cc
# foreach( file_i ${my_srcs})
#     add_custom_command(
#     OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${file_i}
#     COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${file_i} ${CMAKE_CURRENT_BINARY_DIR}/${file_i}
# )
# endforeach( file_i )
       
add_custom_target(wcs_actorgen
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/actor.gen.cc)
  # DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/actor.gen.cc ${CMAKE_CURRENT_BINARY_DIR}/${my_srcs})
  
# add_dependencies(wcs wcs_actorgen)
# Add the header and source files for this directory
set_full_path(THIS_DIR_HEADERS
  actor.hh
  Checkpoint.hh
  Funnel.hh
  simtime.hh
  )

set_full_path(THIS_DIR_SOURCES 
  actor.cc
  )

# # Add the subdirectories

set_source_files_properties(
  ${CMAKE_CURRENT_BINARY_DIR}/actor.gen.cc
  PROPERTIES GENERATED TRUE
)

# set(THIS_DIR_SOURCES "${THIS_DIR_SOURCES}" "${CMAKE_CURRENT_BINARY_DIR}/actor.gen.cc")  

# # Propagate the files up the tree
set(WCS_HEADERS "${WCS_HEADERS}" "${THIS_DIR_HEADERS}" PARENT_SCOPE)
set(WCS_SOURCES "${WCS_SOURCES}" "${THIS_DIR_SOURCES}" PARENT_SCOPE)


