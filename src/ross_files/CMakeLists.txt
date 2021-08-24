# blt_add_executable(
#     NAME ross_files
#     SOURCES actor.cc
#     DEPENDS_ON ROSS mpi
# )
# find_package(Python COMPONENTS Interpreter)

ADD_CUSTOM_COMMAND(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/actor.gen.cc
                   COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/generateInfrastructure.py ${actor_srcs} > ${CMAKE_CURRENT_BINARY_DIR}/actor.gen.cc
                   MAIN_DEPENDENCY generateInfrastructure.py
                   DEPENDS generateInfrastructure.py ${actor_srcs}
                   VERBATIM)
SET_PROPERTY(SOURCE ${CMAKE_CURRENT_BINARY_DIR}/actor.gen.cc
             PROPERTY INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR})



# Add the header and source files for this directory
# set_full_path(THIS_DIR_HEADERS
#   actor.hh
#   Checkpoint.hh
#   Funnel.hh
#   simtime.hh
#   )

# set_full_path(THIS_DIR_SOURCES
#   actor.cc
#   )

# # Add the subdirectories

# # Propagate the files up the tree
# set(WCS_HEADERS "${WCS_HEADERS}" "${THIS_DIR_HEADERS}" PARENT_SCOPE)
# set(WCS_SOURCES "${WCS_SOURCES}" "${THIS_DIR_SOURCES}" PARENT_SCOPE)


