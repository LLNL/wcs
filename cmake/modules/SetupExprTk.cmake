# Download the ExprTk library which consists of a single header file only.
# Set up the variables EXPRTK_DIR and EXPRTK_HEADER.
# Create a target for download step EXPRTK-download

set(WCS_HAS_EXPRTK TRUE)
set(EXPRTK_SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/exprtk)
set(EXPRTK_HEADER exprtk.hpp)

find_file(EXPRTK ${EXPRTK_HEADER}
          HINTS ${EXPRTK_SOURCE_DIR}
                $ENV{EXPRTK_ROOT} ${EXPRTK_ROOT}
                $ENV{EXPRTK_DIR} ${EXPRTK_DIR})

if (EXPRTK)
  message(STATUS "Found ExprTk: ${EXPRTK}")
  add_custom_target(EXPRTK-download)
  unset(EXPRTK_DIR CACHE)
  get_filename_component(EXPRTK_DIR ${EXPRTK} DIRECTORY CACHE)
else ()
  message(STATUS "ExprTk will be downloaded.")
  ExternalProject_Add(EXPRTK
    GIT_REPOSITORY https://github.com/ArashPartow/exprtk.git
    SOURCE_DIR "${EXPRTK_SOURCE_DIR}"
    LOG_DOWNLOAD ON
    STEP_TARGETS download
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )

  unset(EXPRTK_DIR CACHE)
  set(EXPRTK_DIR ${EXPRTK_SOURCE_DIR})

  ExternalProject_Add_StepDependencies(EXPRTK build EXPRTK-download)
endif ()
