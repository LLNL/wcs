# Download the IHeap library which consists of a single header file only.
# Set up the variables IHEAP_DIR and IHEAP_HEADER.
# Create a target for download step IHEAP-download

set(WCS_HAS_IHEAP TRUE)
set(IHEAP_SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/iheap)
set(IHEAP_HEADER iheap.hpp)

find_file(IHEAP ${IHEAP_HEADER}
          HINTS ${IHEAP_SOURCE_DIR}
                $ENV{IHEAP_ROOT} ${IHEAP_ROOT})

if (IHEAP)
  message(STATUS "Found iheap: ${IHEAP}")
  add_custom_target(IHEAP-download)
  unset(IHEAP_DIR CACHE)
  get_filename_component(IHEAP_DIR ${IHEAP} DIRECTORY CACHE)
else ()
  message(STATUS "iheap will be downloaded.")
  ExternalProject_Add(IHEAP
    GIT_REPOSITORY https://github.com/yangle/iheap.git
    SOURCE_DIR "${IHEAP_SOURCE_DIR}"
    LOG_DOWNLOAD ON
    STEP_TARGETS download
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_PATCH ON
    #PATCH_COMMAND "patch ${IHEAP_SOURCE_DIR}/iheap.h < ${CMAKE_SOURCE_DIR}/external/patch_iheap_h.txt"
  )

  unset(IHEAP_DIR CACHE)
  set(IHEAP_DIR ${IHEAP_SOURCE_DIR})

  ExternalProject_Add_StepDependencies(IHEAP build IHEAP-download)
endif ()
