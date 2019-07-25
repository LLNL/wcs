# Download the Catch2 library which consists of a single header file only.
# Set up the variables CATCH2_DIR and CATCH2_HEADER.
# Create a target for download step CATCH2-download

set(WCS_HAS_CATCH2 TRUE)
set(CATCH2_SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/Catch2)
set(CATCH2_HEADER catch.hpp)

find_file(CATCH2 ${CATCH2_HEADER}
          HINTS ${CATCH2_SOURCE_DIR}
                $ENV{CATCH2_ROOT} ${CATCH2_ROOT})

if (CATCH2)
  message(STATUS "Found Catch2: ${CATCH2}")
  add_custom_target(CATCH2-download)
  unset(CATCH2_DIR CACHE)
  get_filename_component(CATCH2_DIR ${CATCH2} DIRECTORY CACHE)
else ()
  message(STATUS "Catch2 not found")
  ExternalProject_Add(CATCH2
    URL http://github.com/catchorg/Catch2/releases/download/v2.9.1/catch.hpp 
    DOWNLOAD_DIR "${CATCH2_SOURCE_DIR}"
    DOWNLOAD_NAME "${CATCH2_HEADER}"
    DOWNLOAD_NO_EXTRACT TRUE
    LOG_DOWNLOAD ON
    STEP_TARGETS download
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )

  unset(CATCH2_DIR CACHE)
  set(CATCH2_DIR ${CATCH2_SOURCE_DIR})

  ExternalProject_Add_StepDependencies(CATCH2 build CATCH2-download)
endif ()
