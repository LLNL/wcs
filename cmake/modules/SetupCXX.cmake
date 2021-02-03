# This sets up all of the proper compiler information, including some
# custom flags. <Tom's skeptical face>

include(CheckCXXCompilerFlag)
include(CheckIncludeFileCXX)

# MACRO WCS_CHECK_AND_APPEND_FLAG
#
# Purpose: checks that all flags are valid and appends them to the
#   given list. Valid means that the compiler does not throw an error
#   upon encountering the flag.
#
# Arguments:
#   VAR The list of current flags
#   ARGN The flags to check
#
# Note: If flag is not valid, it is not appended.
macro(wcs_check_and_append_flag MY_VAR)
  foreach(flag ${ARGN})
    string(REPLACE "-" "_" _CLEAN_FLAG "${flag}")

    set(CMAKE_REQUIRED_LIBRARIES "${flag}")
    check_cxx_compiler_flag("${flag}" FLAG_${_CLEAN_FLAG}_OK)
    unset(CMAKE_REQUIRED_LIBRARIES)

    if (FLAG_${_CLEAN_FLAG}_OK)
      set(${MY_VAR} "${${MY_VAR}} ${flag}")
    endif ()
  endforeach()
endmacro()

# Temporary workaround to force CMake to recognize the XL
# compiler. The latest beta XL compiler is recognized as Clang.
if (CMAKE_CXX_COMPILER MATCHES "xlc" AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_COMPILER_ID "XL")
endif ()


################################################################
# Muck with flags
################################################################

if (${WCS_GPROF})
  set (COMPILER_OPT_FOR_GPROF "-pg")
  set (WCS_PERF_PROF ON)
endif (${WCS_GPROF})

# Initialize C++ flags
wcs_check_and_append_flag(CMAKE_CXX_FLAGS
  -fPIC -g -Wall -Wextra -Wno-unused-parameter -Wnon-virtual-dtor
  -Wno-deprecated-declarations -std=c++17 ${COMPILER_OPT_FOR_GPROF} ${USER_FLAGS})
  #taking out -Wshadow as ExprTK generates too much warnings

# Disable all optimization in debug for better viewing under debuggers
# (cmake already adds -g)
wcs_check_and_append_flag(CMAKE_CXX_FLAGS_DEBUG -O0)

if (${UPPER_PROJECT_NAME}_WARNINGS_AS_ERRORS)
  wcs_check_and_append_flag(_WERROR_FLAGS -Werror)
  separate_arguments(_WERROR_FLAGS NATIVE_COMMAND "${_WERROR_FLAGS}")
  if (NOT TARGET WCS_CXX_FLAGS_werror)
    add_library(WCS_CXX_FLAGS_werror INTERFACE)
    set_property(TARGET WCS_CXX_FLAGS_werror PROPERTY
      INTERFACE_COMPILE_OPTIONS $<$<COMPILE_LANGUAGE:CXX>:${_WERROR_FLAGS}>)

    # Add the "library" to the export
    install(TARGETS WCS_CXX_FLAGS_werror EXPORT WCSTargets)
  endif ()
endif ()

# Some behavior is dependent on the compiler version.
if (NOT CMAKE_CXX_COMPILER_VERSION)
  execute_process(
    COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
    OUTPUT_VARIABLE CXX_VERSION)
else ()
  set(CXX_VERSION "${CMAKE_CXX_COMPILER_VERSION}")
endif ()

# - Special handling if we're compiling with Clang's address sanitizer
# - gcc toolchain handling for interoperability, especially with the exteral
#   libraries pre-built using gcc
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  if (USE_CLANG_LIBCXX)
    wcs_check_and_append_flag(CMAKE_CXX_FLAGS "--stdlib=libc++")
  else (USE_CLANG_LIBCXX)
    if (USE_GCC_LIBCXX)
      wcs_check_and_append_flag(CMAKE_CXX_FLAGS
        "--gcc-toolchain=/usr/tce/packages/gcc/gcc-${GCC_TOOLCHAIN_VER}")
    endif (USE_GCC_LIBCXX)
  endif (USE_CLANG_LIBCXX)

  if (CMAKE_BUILD_TYPE MATCHES Debug)
    wcs_check_and_append_flag(CMAKE_CXX_FLAGS
      -fsanitize=address -fno-omit-frame-pointer -fsanitize-recover=address)
  else()
    wcs_check_and_append_flag(CMAKE_CXX_FLAGS -fno-omit-frame-pointer)
  endif ()
endif ()

# Turn off some annoying warnings
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  # Bugs with Intel compiler version 19
  #https://community.intel.com/t5/Intel-C-Compiler/quot-if-constexpr-quot-and-quot-missing-return-statement-quot-in/td-p/1154551
  #https://bitbucket.org/berkeleylab/upcxx/issues/286/icc-bug-bogus-warning-use-of-offsetof-with
  # set(GCC_PATH "-gcc-name=/usr/tce/packages/gcc/gcc-8.1.0/bin/gcc")
  if (GCC_PATH)
    set(GCC_INTEROP "-gcc-name=${GCC_PATH}")
  endif (GCC_PATH)
  wcs_check_and_append_flag(CMAKE_CXX_FLAGS -diag-disable=2196 -wd1011 -wd1875 ${GCC_INTEROP})

endif ()


################################################################
# Check if std::filesystem is available
################################################################
try_compile(WCS_HAS_STD_FILESYSTEM "${CMAKE_BINARY_DIR}/temp"
            "${CMAKE_SOURCE_DIR}/cmake/tests/has_filesystem.cpp"
            CMAKE_FLAGS ${CMAKE_CXX_FLAGS}
            LINK_LIBRARIES stdc++fs)
if (WCS_HAS_STD_FILESYSTEM)
  message(STATUS "Compiler has std::filesystem support")
else ()
  message(STATUS "Compiler does not have std::filesystem support. Using boost::filesystem")
endif (WCS_HAS_STD_FILESYSTEM)


################################################################
# Initialize RPATH (always full RPATH)
# Note: see https://cmake.org/Wiki/CMake_RPATH_handling
################################################################

# Use RPATH on OS X
if (APPLE)
  set(CMAKE_MACOSX_RPATH ON)
endif ()

# Use (i.e. don't skip) RPATH for build
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# Use same RPATH for build and install
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# Add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
  "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" _IS_SYSTEM_DIR)
if (${_IS_SYSTEM_DIR} STREQUAL "-1")
  # Set the install RPATH correctly
  list(APPEND CMAKE_INSTALL_RPATH
    "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif ()

set(CMAKE_C_COMPILER_ID ${CMAKE_CXX_COMPILER_ID})

if (WCS_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
endif (WCS_WITH_OPENMP)

# Testing for compiler feature supports 
#include(CheckCXXSourceCompiles)
