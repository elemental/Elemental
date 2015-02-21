#
#  Copyright 2009-2015, Jack Poulson
#  All rights reserved.
#
#  Copyright 2012, Jed Brown
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
cmake_minimum_required(VERSION 2.8.8)
include(CheckCSourceCompiles)
include(CheckCXXSourceCompiles)
include(CheckFunctionExists)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/external/cmake/")
include(GetGitRevisionDescription)
include(language_support_v2)

project(Elemental)
set(EL_VERSION_MAJOR 0)
set(EL_VERSION_MINOR 86-dev)

# Configuration options
# =====================

# Basic options
# -------------

# "CMAKE_BUILD_TYPE"
#   Elemental requires it to be one of the following choices:
#     1. "PureDebug": Vanilla MPI build meant for development and debugging
#     2. "PureRelease": Vanilla MPI build meant for production runs
#     3. "HybridDebug": MPI+OpenMP build meant for development and debugging
#     4. "HybridRelease": MPI+OpenMP build meant for production runs
#   If "CMAKE_BUILD_TYPE" is not defined, the default is "PureRelease".

# "CXX_FLAGS"
#   Optimization and debugging/symbol flags

# "OpenMP_CXX_FLAGS"
#   Overrides the default compile flags for adding OpenMP support to CXX code

# TODO: Add list of MPI variables here, such as "MPI_CXX_COMPILE_FLAGS"

# TODO: Add discussion of "MATH_LIBS" and "REFERENCE_ROOT"

set(BUILD_SHARED_LIBS "Build shared libraries?" ON)

option(EL_C_INTERFACE "Build C interface" ON)
if(BUILD_SHARED_LIBS AND EL_C_INTERFACE)
  # Define PYTHON_SITE_PACKAGES if you want to install the python package
  # somewhere other than the default system-wide location (e.g., within your
  # home directory)
  option(INSTALL_PYTHON_PACKAGE "Install python interface?" TRUE)
endif()

# Whether or not to use Qt5 for visualization if it is found.
# Setting Qt5Widgets_DIR, Qt5_DIR, or Qt5_LIBDIR overrides this variable
# and is often necessary anyway.
option(EL_USE_QT5 "Attempt to use Qt5?" OFF)

# NOTE: This is not yet used within Elemental...
#
# If BUILD_KISSFFT is set to 'OFF', then an external KISSFFT can be
# built by defining:
#   KISSFFT_INCLUDE_DIR: The directory containing kiss_fft.h
#   KISSFFT_LIBS: The list of KISS-FFT libraries to link
option(BUILD_KISSFFT "Build KISS FFT" OFF)

option(EL_EXAMPLES "Build simple examples?" OFF)
option(EL_TESTS "Build performance and correctness tests?" OFF)
option(EL_EXPERIMENTAL "Build experimental code" OFF)

# Whether or not to have the Memory class zero initialize what it allocates.
# If valgrind was detected and is running, this will be forced anyway.
option(EL_ZERO_INIT "Initialize buffers to zero by default?" OFF)

# Attempt to use 64-bit integers
option(EL_USE_64BIT_INTS "Use 64-bit integers where possible" OFF)

option(EL_USE_CUSTOM_ALLTOALLV "Avoid MPI_Alltoallv for performance reasons" ON)
option(EL_BARRIER_IN_ALLTOALLV "Barrier before posting non-blocking recvs" OFF)

# If the version of METIS packaged with Elemental is to be built (the default),
# then no METIS-specific variables need to be specified, but if the user prefers
# to use their own version, then the root path of the installation should be
# passed in via the METIS_ROOT variable.
option(BUILD_METIS "Build METIS" ON)

# ParMETIS cannot yet be legally distributed under an unrestrictive 
# license and is currently limited to educational and research purposes.
# (recent versions of METIS can be distributed under the Apache License 
# Version 2).
#
# However, if the user is in an educational or government institution, then 
# it is possible to re-enable support for ParMETIS by specifying the following
# variables:
#
#   METIS_ROOT: the full path to the root of the MeTiS source tree
#   PARMETIS_ROOT: the full path to the root of the ParMeTiS source tree
#   PARMETIS_LIBS: the list of ParMETIS libraries to link to
#   PARMETIS_TLS_PATCH: defined iff there is a GKlib configure file
#                       (if you have not heard of this, ignore it)
#   GKLIB_INCLUDE_DIR: if PARMETIS_TLS_PATCH is defined, this should point
#                      to the directory containing gklib_tls.h
# 
# and adding in a custom ParMETIS routine for performing nodal bisections
# in parallel (which was originally part of Clique and can be easily found).
# Please contact dev@libelemental.org or poulson@stanford.edu for more details.
#
#option(BUILD_PARMETIS "Build ParMETIS" OFF)

# Advanced options
# ----------------

# Since it is surprisingly common for MPI libraries to have bugs in their 
# support for complex data, the following option forces Elemental to cast 
# all possible MPI communications in terms of twice as many real units of data.
option(EL_AVOID_COMPLEX_MPI "Avoid potentially buggy complex MPI routines" ON)
mark_as_advanced(EL_AVOID_COMPLEX_MPI)

# At one point, a bug was found in IBM's C++ compiler for Blue Gene/P, 
# where OpenMP statements of the form a[i] += alpha b[i], with complex data,
# would segfault and/or return incorrect results
option(EL_AVOID_OMP_FMA "Avoid a bug in the IBM compilers." OFF)
mark_as_advanced(EL_AVOID_OMP_FMA)

# Due to a subtle flaw in the Blue Gene/P extensions for MPICH2, treating 
# floating-point data as a collection of byte-sized objects results in a 
# better algorithm being chosen for MPI_Allgather. This should not effect
# performance on most machines.
option(EL_USE_BYTE_ALLGATHERS "Avoid BG/P allgather performance bug." ON)
mark_as_advanced(EL_USE_BYTE_ALLGATHERS)

# If MPI_Reduce_scatter_block doesn't exist, perform it by composing 
# MPI_Allreduce and std::memcpy rather than MPI_Reduce and MPI_Scatter
option(EL_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
       "AllReduce based block MPI_Reduce_scatter" OFF)
mark_as_advanced(EL_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE)

# Print a warning any time a redistribution is performed which unpacks a 
# large amount of data with a non-unit stride
option(EL_CACHE_WARNINGS "Warns when using cache-unfriendly routines" OFF)
mark_as_advanced(EL_CACHE_WARNINGS)

# Print a warning when an improperly aligned redistribution is performed, 
# i.e., if an unnecessary permutation communication stage must take place
option(EL_UNALIGNED_WARNINGS 
       "Warn when performing unaligned redistributions" OFF)
mark_as_advanced(EL_UNALIGNED_WARNINGS)

# Print a warning if an opportunity was missed to implement a redistribution
# approach specifically for vectors (instead of matrices)
option(EL_VECTOR_WARNINGS 
       "Warn when vector redistribution chances are missed" OFF)
mark_as_advanced(EL_VECTOR_WARNINGS)

# Handle RPATHs for Mac
# =====================
# These settings are recommend by 
# <http://www.cmake.org/Wiki/CMake_RPATH_handling>, which was accessed on
# February 20, 2015.

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

# Build logic
# ===========

# Elemental must be built "out-of-source", so we start by ensuring that the
# source and build directories are different.
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "In-source build attempted; please clean the CMake cache and then switch to an out-of-source build, e.g., rm CMakeCache.txt && rm -Rf CMakeFiles/ && mkdir build/ && cd build/ && cmake ..")
endif()

# Get the Git revision
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

# Extract a few booleans from the build type and default to PureRelease
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Build mode not specified, defaulting to PureRelease build.")
  set(CMAKE_BUILD_TYPE PureRelease)
endif()
if(CMAKE_BUILD_TYPE STREQUAL "HybridDebug")
  set(EL_HYBRID TRUE)
  set(EL_DEBUG TRUE)
elseif(CMAKE_BUILD_TYPE STREQUAL "HybridRelease")
  set(EL_HYBRID TRUE)
  set(EL_RELEASE TRUE)
elseif(CMAKE_BUILD_TYPE STREQUAL "PureDebug")
  set(EL_PURE TRUE)
  set(EL_DEBUG TRUE)
elseif(CMAKE_BUILD_TYPE STREQUAL "PureRelease")
  set(EL_PURE TRUE)
  set(EL_RELEASE TRUE)
else()
  message(FATAL_ERROR "CMAKE_BUILD_TYPE must be [Hybrid,Pure][Debug,Release]")
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} UPPER_BUILD_TYPE)

# Set the basic compile flags from the build type
if(NOT WIN32) 
  set(LANGUAGES CXX C Fortran)
  foreach(LANG ${LANGUAGES})
    if(NOT ${LANG}_FLAGS)
      if(EL_DEBUG)
        #set(${LANG}_FLAGS "-O2 -g -Wunused-variable -Wunused-but-set-variable -Wunused-local-typedefs" CACHE STRING
        set(${LANG}_FLAGS "-O2 -g" CACHE STRING
          "${LANG} optimization/debug flags for ${UPPER_BUILD_TYPE} mode")
      else()
        set(${LANG}_FLAGS "-O3" CACHE STRING
          "${LANG} optimization flags for ${UPPER_BUILD_TYPE} mode")
      endif()
    endif()
  endforeach()
endif()

if(BUILD_PARMETIS)
  message(FATAL_ERROR "ParMETIS cannot be legally distributed with Elemental under a permissive license. A version of Elemental restricted to educational purposes which contains ParMETIS will hopefully soon be available.")
#  add_subdirectory(external/parmetis)
#  set(EL_HAVE_METIS TRUE)
#  set(EL_HAVE_PARMETIS TRUE)
#
#  set(PARMETIS_ROOT ${PROJECT_SOURCE_DIR}/external/parmetis)
#  set(METIS_ROOT ${PROJECT_SOURCE_DIR}/external/parmetis/metis)
#  set(GKLIB_INCLUDE_DIR ${PROJECT_BINARY_DIR}/external/parmetis/include)
#
#  include_directories(${GKLIB_INCLUDE_DIR})
#  include_directories(${METIS_ROOT}/include)
#  include_directories(${PARMETIS_ROOT}/include)
elseif(METIS_ROOT AND PARMETIS_ROOT)
  message(FATAL_ERROR "ParMETIS cannot be legally distributed with Elemental under a permissive license. A version of Elemental restricted to educational purposes which contains ParMETIS will hopefully soon be available.")
#  set(EL_HAVE_METIS TRUE)
#  set(EL_HAVE_PARMETIS TRUE)
#  if(PARMETIS_TLS_PATCH)
#    include_directories(${GKLIB_INCLUDE_DIR})
#  endif()
#
#  include_directories(${METIS_ROOT}/include)
#  include_directories(${PARMETIS_ROOT}/include)
#  set(MANUAL_PARMETIS TRUE)
elseif(BUILD_METIS)
  add_subdirectory(external/metis)

  set(METIS_ROOT ${PROJECT_SOURCE_DIR}/external/metis)
  include_directories(${METIS_ROOT}/include)

  set(EL_HAVE_METIS TRUE)
else()
  if(NOT METIS_ROOT)
    message(FATAL_ERROR 
      "If METIS is not to be built, then METIS_ROOT must be specified")
  endif()
  include_directories(${METIS_ROOT}/include)

  set(EL_HAVE_METIS TRUE)
  set(MANUAL_METIS TRUE)
endif()

if(BUILD_KISSFFT)
  add_subdirectory(external/kiss_fft)
  include_directories(external/kiss_fft)
elseif(KISSFFT_INCLUDE_DIR AND KISSFFT_LIBS)
  include_directories(${KISSFFT_INCLUDE_DIR})
  set(MANUAL_KISSFFT TRUE)
endif()

if(BUILD_SHARED_LIBS AND WIN32)
  option(EL_DECLSPEC "Use __declspec(dllexport)?" ON)
else()
  option(EL_DECLSPEC "Use __declspec(dllexport)?" OFF)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(tests/Fortran)
include(tests/CXX)
include(tests/MPI)
include(tests/OpenMP)

include(tests/Qt5)
if(EL_HAVE_QT5)
  set(EL_HEADERS_PREMOC 
      "include/El/io/DisplayWindow-premoc.hpp;include/El/io/ComplexDisplayWindow-premoc.hpp")
  qt_wrap_cpp(El EL_MOC_SRC ${EL_HEADERS_PREMOC})
  include_directories(${Qt5Widgets_INCLUDE_DIRS})
  add_definitions(${Qt5Widgets_DEFINITIONS})
  set(EXTRA_FLAGS "${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS} ${EXTRA_FLAGS}")

  # Qt5Widgets_DIR = Qt5_LIBDIR/cmake/Qt5Widgets
  get_filename_component(Qt5_CMAKEDIR ${Qt5Widgets_DIR} PATH)
  get_filename_component(Qt5_LIBDIR ${Qt5_CMAKEDIR} PATH)
endif()

include(FindValgrind)
if(VALGRIND_FOUND)
  include_directories(${VALGRIND_INCLUDE_DIR})
  set(EL_HAVE_VALGRIND TRUE)
endif()

include(tests/Math)

# Now append the extra flags
set(CXX_FLAGS "${CXX_FLAGS} ${EXTRA_FLAGS}")
set(CMAKE_CXX_FLAGS_${UPPER_BUILD_TYPE} ${CXX_FLAGS})

# Add the Parallel Multiple Relatively Robust Representations (PMRRR) project
add_subdirectory(external/pmrrr)

# Create the Elemental configuration header
configure_file(${PROJECT_SOURCE_DIR}/cmake/config.h.cmake
               ${PROJECT_BINARY_DIR}/include/El/config.h)
install(FILES ${PROJECT_BINARY_DIR}/include/El/config.h 
        DESTINATION include/El)

# Create a file which can be included in Makefile's.
# This is meant to be analogous to PETSc's 'conf/petscvariables' file
set(MPI_CXX_INCSTRING)
foreach(INC_PATH ${MPI_CXX_INCLUDE_PATH})
  set(MPI_CXX_INCSTRING "${MPI_CXX_INCSTRING} -I${INC_PATH}")  
endforeach()
set(MATH_LIBSTRING)
foreach(LIB ${MATH_LIBS})
  set(MATH_LIBSTRING "${MATH_LIBSTRING} ${LIB}")
endforeach()
set(MPI_CXX_LIBSTRING)
foreach(LIB ${MPI_CXX_LIBRARIES})
  set(MPI_CXX_LIBSTRING "${MPI_CXX_LIBSTRING} ${LIB}")
endforeach()
set(METIS_LIBSTRING)
foreach(LIB ${METIS_LIBS})
  set(METIS_LIBSTRING "${METIS_LIBSTRING} ${LIB}")
endforeach()
set(QT5_DEFSTRING)
foreach(DEF ${Qt5Widgets_DEFINITIONS})
  set(QT5_DEFSTRING "${QT5_DEFSTRING} ${DEF}")
endforeach()
set(QT5_INCSTRING)
foreach(INC ${Qt5Widgets_INCLUDE_DIRS})
  set(QT5_INCSTRING "${QT5_INCSTRING} -I${INC}")
endforeach()
set(QT5_COMPILESTRING)
foreach(FLAG ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS})
  set(QT5_COMPILESTRING "${QT5_COMPILE_STRING} ${FLAG}")
endforeach()
# TODO: Generalize this for non-Unix architectures
set(QT5_LIBSTRING "-L${Qt5_LIBDIR} -lQt5Widgets -lQt5Gui -lQt5Core")
configure_file(${PROJECT_SOURCE_DIR}/cmake/ElVars.cmake
               ${PROJECT_BINARY_DIR}/conf/ElVars @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/conf/ElVars DESTINATION conf)

# The main library
# ================

# Grab the .c and .cpp files
file(GLOB_RECURSE EL_C_CPP_SOURCE RELATIVE ${PROJECT_SOURCE_DIR} 
  "src/*.c" "src/*.cpp")
# Grab the C/C++ headers
file(GLOB_RECURSE EL_C_CPP_HEADERS RELATIVE ${PROJECT_SOURCE_DIR} 
  "include/*.h" "include/*.hpp")
set(EL_C_CPP_FILES "${EL_C_CPP_SOURCE};${EL_C_CPP_HEADERS};${EL_MOC_SRC}")

set(LINK_LIBS ${MATH_LIBS} ${MPI_CXX_LIBRARIES})

if(EL_HAVE_QT5)
  set(LINK_LIBS ${LINK_LIBS} ${Qt5Widgets_LIBRARIES})  
endif()

if(BUILD_KISSFFT)
  set(LINK_LIBS ${LINK_LIBS} kiss_fft)
elseif(MANUAL_KISSFFT)
  set(LINK_LIBS ${LINK_LIBS} ${KISSFFT_LIBS})
endif()

if(BUILD_PARMETIS)
  set(LINK_LIBS ${LINK_LIBS} parmetis)
elseif(MANUAL_PARMETIS)
  set(LINK_LIBS ${LINK_LIBS} ${PARMETIS_LIBS})
elseif(BUILD_METIS)
  set(LINK_LIBS ${LINK_LIBS} metis)
elseif(MANUAL_METIS)
  set(LINK_LIBS ${LINK_LIBS} ${METIS_LIBS})
endif()

add_library(El ${EL_C_CPP_FILES})
target_link_libraries(El pmrrr ${LINK_LIBS})
if(MPI_LINK_FLAGS)
  set_target_properties(El PROPERTIES LINK_FLAGS ${MPI_LINK_FLAGS})
endif()
install(TARGETS El DESTINATION lib)

# Define the header-file preparation rules
set(PREPARED_HEADERS)
foreach(HEADER ${EL_C_CPP_HEADERS})
  add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/${HEADER}
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/${HEADER}
            ${PROJECT_BINARY_DIR}/${HEADER}
    DEPENDS "${PROJECT_SOURCE_DIR}/${HEADER}")
  list(APPEND PREPARED_HEADERS ${PROJECT_BINARY_DIR}/${HEADER})

  get_filename_component(HEADER_PATH ${HEADER} PATH)
  install(FILES ${PROJECT_BINARY_DIR}/${HEADER} DESTINATION ${HEADER_PATH})
endforeach()
add_custom_target(prepare_El_headers DEPENDS ${PREPARED_HEADERS})
add_dependencies(El prepare_El_headers)

if(BUILD_SHARED_LIBS)
  if(INSTALL_PYTHON_PACKAGE)
    if(NOT PYTHON_SITE_PACKAGES)
      execute_process( 
        COMMAND python -c 
        "from distutils.sysconfig import get_python_lib; print get_python_lib()" 
        OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    endif()
    message(STATUS "PYTHON_SITE_PACKAGES=${PYTHON_SITE_PACKAGES}")
    install(
      DIRECTORY ${PROJECT_SOURCE_DIR}/include/
      DESTINATION ${PYTHON_SITE_PACKAGES}
      FILES_MATCHING PATTERN "*.py"
    )
  else()
    install(
      DIRECTORY ${PROJECT_SOURCE_DIR}/include/
      DESTINATION ${CMAKE_INSTALL_PREFIX}/include/
      FILES_MATCHING PATTERN "*.py"
    )
  endif()
endif()

# TODO: Ensure that source files depend on header files

# Make sure the Elemental headers can be found
include_directories("${PROJECT_BINARY_DIR}/include")

# Build the test drivers if necessary
if(EL_TESTS)
  set(TEST_DIR ${PROJECT_SOURCE_DIR}/tests)
  set(TEST_TYPES core blas_like lapack_like optimization)
  foreach(TYPE ${TEST_TYPES})
    file(GLOB_RECURSE ${TYPE}_TESTS 
      RELATIVE ${PROJECT_SOURCE_DIR}/tests/${TYPE}/ "tests/${TYPE}/*.cpp")

    set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/tests/${TYPE}")
    foreach(TEST ${${TYPE}_TESTS})
      set(DRIVER ${TEST_DIR}/${TYPE}/${TEST})
      get_filename_component(TESTNAME ${TEST} NAME_WE)
      add_executable(tests-${TYPE}-${TESTNAME} ${DRIVER})
      set_source_files_properties(${DRIVER} PROPERTIES 
        OBJECT_DEPENDS "${PREPARED_HEADERS}")
      target_link_libraries(tests-${TYPE}-${TESTNAME} El)
      set_target_properties(tests-${TYPE}-${TESTNAME} PROPERTIES
        OUTPUT_NAME ${TESTNAME} RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
      if(MPI_LINK_FLAGS)
        set_target_properties(tests-${TYPE}-${TESTNAME} PROPERTIES
          LINK_FLAGS ${MPI_LINK_FLAGS})
      endif()
      install(TARGETS tests-${TYPE}-${TESTNAME} DESTINATION bin/tests/${TYPE})
    endforeach()
  endforeach()
endif()

# Build the example drivers if necessary
if(EL_EXAMPLES)
  set(EXAMPLE_DIR ${PROJECT_SOURCE_DIR}/examples)
  set(EXAMPLE_TYPES 
      optimization core blas_like interface io lapack_like matrices)
  foreach(TYPE ${EXAMPLE_TYPES})
    if(EL_C_INTERFACE)
      file(GLOB_RECURSE ${TYPE}_EXAMPLES RELATIVE 
           ${PROJECT_SOURCE_DIR}/examples/${TYPE}/ "examples/${TYPE}/*.cpp"
                                                   "examples/${TYPE}/*.c" )
    else()
      file(GLOB_RECURSE ${TYPE}_EXAMPLES RELATIVE 
           ${PROJECT_SOURCE_DIR}/examples/${TYPE}/ "examples/${TYPE}/*.cpp")
    endif()
    set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/examples/${TYPE}")
    foreach(EXAMPLE ${${TYPE}_EXAMPLES})
      set(DRIVER ${EXAMPLE_DIR}/${TYPE}/${EXAMPLE})
      get_filename_component(EXNAME ${EXAMPLE} NAME_WE)
      add_executable(examples-${TYPE}-${EXNAME} ${DRIVER})
      set_source_files_properties(${DRIVER} PROPERTIES 
        OBJECT_DEPENDS "${PREPARED_HEADERS}")
      target_link_libraries(examples-${TYPE}-${EXNAME} El)
      set_target_properties(examples-${TYPE}-${EXNAME} PROPERTIES 
        OUTPUT_NAME ${EXNAME} RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
      if(MPI_LINK_FLAGS)
        set_target_properties(examples-${TYPE}-${EXNAME} PROPERTIES
          LINK_FLAGS ${MPI_LINK_FLAGS})
      endif()
      install(TARGETS examples-${TYPE}-${EXNAME} 
        DESTINATION bin/examples/${TYPE})
    endforeach()
  endforeach()
endif()

# Build experimental drivers
if(EL_EXPERIMENTAL)
  set(EXPERIMENTAL_DIR ${PROJECT_SOURCE_DIR}/experimental)

  # Build the G3D example(s)
  set(G3D_EXPERS G3DGemm)
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/experimental/g3d")
  foreach(EXPER ${G3D_EXPERS})
    set(DRIVER ${EXPERIMENTAL_DIR}/g3d/${EXPER}.cpp)
    add_executable(experimental-g3d-${EXPER} ${DRIVER})
    set_source_files_properties(${DRIVER} PROPERTIES 
      OBJECT_DEPENDS "${PREPARED_HEADERS}")
    target_link_libraries(experimental-g3d-${EXPER} El)
    set_target_properties(experimental-g3d-${EXPER} PROPERTIES 
      OUTPUT_NAME ${EXPER} RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(experimental-g3d-${EXPER} PROPERTIES
        LINK_FLAGS ${MPI_LINK_FLAGS})
    endif()
    install(TARGETS experimental-g3d-${EXPER} DESTINATION bin/experimental/g3d)
  endforeach()
endif()

# If Elemental is a subproject, then pass some variables to the parent
if(NOT CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  include(./cmake/ElSub.cmake)
endif()
