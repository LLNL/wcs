# Download the Catch2 library which consists of a single header file only.
# Set up the variables CATCH2_INCLUDE_DIR and CATCH2_HEADER.
# Create a target for download step CATCH2

set(WCS_HAS_CATCH2 TRUE)
set(CATCH2_SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/Catch2)
set(CATCH2_HEADER catch.hpp)

find_file(CATCH2 ${CATCH2_HEADER}
          HINTS ${CATCH2_SOURCE_DIR}
                $ENV{CATCH2_ROOT} ${CATCH2_ROOT}
          PATH_SUFFIXES single_include/catch2)

if (CATCH2)
  message(STATUS "Found Catch2: ${CATCH2}")
  add_custom_target(CATCH2)
  get_filename_component(CATCH2_DIR ${CATCH2} DIRECTORY)
  get_filename_component(CATCH2_INCLUDE_DIR ${CATCH2_DIR} DIRECTORY CACHE)
else ()
  message(STATUS "Catch2 not found")
  ExternalProject_Add(CATCH2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    SOURCE_DIR ${CATCH2_SOURCE_DIR}
    LOG_DOWNLOAD ON
    CMAKE_ARGS -DCATCH_BUILD_TESTING:BOOL=OFF -DCATCH_INSTALL_DOCS:BOOL=OFF -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/catch2
  )

  set(CATCH2_INCLUDE_DIR ${CATCH2_SOURCE_DIR}/single_include)
  set(CATCH2_CONTRIB_DIR ${CATCH2_SOURCE_DIR}/contrib)
endif ()
