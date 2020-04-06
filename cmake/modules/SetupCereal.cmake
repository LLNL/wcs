# Download the Cereal library which consists of a single header file only.
# Set up the variables CEREAL_DIR and CEREAL_HEADER.
# Create a target for download step CEREAL-download

set(WCS_HAS_CEREAL TRUE)
set(CEREAL_SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/cereal)
set(CEREAL_HEADER cereal.hpp)

find_file(CEREAL ${CEREAL_HEADER}
          HINTS ${CEREAL_SOURCE_DIR}/include/cereal
                $ENV{CEREAL_ROOT}/include/cereal
                ${CEREAL_ROOT}/include/cereal)

if (CEREAL)
  message(STATUS "Found Cereal: ${CEREAL}")
  add_custom_target(CEREAL-download)
  unset(CEREAL_DIR CACHE)
  get_filename_component(CEREAL_CEREAL_DIR ${CEREAL} DIRECTORY)
  get_filename_component(CEREAL_HEADER_DIR ${CEREAL_CEREAL_DIR} PATH CACHE)
  get_filename_component(CEREAL_DIR ${CEREAL_HEADER_DIR} PATH CACHE)
#  message(STATUS "CEREAL_DIR: ${CEREAL_DIR}")
  message(STATUS "CEREAL_HEADER_DIR: ${CEREAL_HEADER_DIR}")
else ()
  message(STATUS "Cereal will be downloaded.")
  ExternalProject_Add(CEREAL
    GIT_REPOSITORY https://github.com/USCiLab/cereal.git
    SOURCE_DIR "${CEREAL_SOURCE_DIR}"
    LOG_DOWNLOAD ON
    STEP_TARGETS download
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )

  unset(CEREAL_DIR CACHE)
  set(CEREAL_DIR ${CEREAL_SOURCE_DIR})
  set(CEREAL_HEADER_DIR ${CEREAL_SOURCE_DIR}/include)
#  message(STATUS "CEREAL_DIR: ${CEREAL_DIR}")
  message(STATUS "CEREAL_HEADER_DIR: ${CEREAL_HEADER_DIR}")

  ExternalProject_Add_StepDependencies(CEREAL build CEREAL-download)
endif ()
