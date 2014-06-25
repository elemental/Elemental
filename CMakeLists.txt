#
#  Copyright (c) 2009-2014, Jack Poulson
#                     2012, Jed Brown 
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
set(Elemental_VERSION_MAJOR 0)
set(Elemental_VERSION_MINOR 84-dev)

################################################################################
# Common build options                                                         #
################################################################################

# By default, Elemental builds static libraries, but shared libraries are 
# also supported through this option.
option(BUILD_SHARED_LIBS "Switch to shared libraries" OFF)

# Whether or not to build a collection of simple example drivers.
option(ELEM_EXAMPLES "A collection of simple examples" OFF)

# Whether or not to build a collection of performance and correctness tests
option(ELEM_TESTS "A collection of performance and correctness tests" OFF)

# Whether or not to force the Memory class zero initialize what it allocates.
# If valgrind was detected and is running, this will be forced anyway.
option(ELEM_ZERO_INIT "Initialize buffers to zero by default?" OFF)

# Whether or not to build current 'experimental' code, such as:
# Martin's 3D Gemm code
option(ELEM_EXPERIMENTAL "Build experimental code" OFF)

# Attempt to use 64-bit integers
option(ELEM_USE_64BIT_INTS "Use 64-bit integers where possible" OFF)

# Whether or not to use Qt5 for visualization if it is found (experimental)
# NOTE: Setting Qt5Widgets_DIR, Qt5_DIR, or Qt5_LIBDIR overrides this variable
option(ELEM_USE_QT5 "Attempt to use Qt5?" OFF)

################################################################################
# Advanced build options                                                       #
################################################################################

# Since it is surprisingly common for MPI libraries to have bugs in their 
# support for complex data, the following option forces Elemental to cast 
# all possible MPI communications in terms of twice as many real units of data.
option(ELEM_AVOID_COMPLEX_MPI "Avoid potentially buggy complex MPI routines" ON)
mark_as_advanced(ELEM_AVOID_COMPLEX_MPI)

# At one point, a bug was found in IBM's C++ compiler for Blue Gene/P, 
# where OpenMP statements of the form a[i] += alpha b[i], with complex data,
# would segfault and/or return incorrect results
option(ELEM_AVOID_OMP_FMA "Avoid a bug in the IBM compilers." OFF)
mark_as_advanced(ELEM_AVOID_OMP_FMA)

# Due to a subtle flaw in the Blue Gene/P extensions for MPICH2, treating 
# floating-point data as a collection of byte-sized objects results in a 
# better algorithm being chosen for MPI_Allgather. This should not effect
# performance on most machines.
option(ELEM_USE_BYTE_ALLGATHERS "Avoid BG/P allgather performance bug." ON)
mark_as_advanced(ELEM_USE_BYTE_ALLGATHERS)

# If MPI_Reduce_scatter_block doesn't exist, perform it by composing 
# MPI_Allreduce and std::memcpy rather than MPI_Reduce and MPI_Scatter
option(ELEM_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
       "AllReduce based block MPI_Reduce_scatter" OFF)
mark_as_advanced(ELEM_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE)

# This can easily be performed, but is likely not a good idea
option(ELEM_POOL_MEMORY 
       "Make Memory class accumulate memory until destruction" OFF)
mark_as_advanced(ELEM_POOL_MEMORY)

################################################################################
# Elemental-development build options                                          #
################################################################################

# Print a warning any time a redistribution is performed which unpacks a 
# large amount of data with a non-unit stride
option(ELEM_CACHE_WARNINGS "Warns when using cache-unfriendly routines" OFF)
mark_as_advanced(ELEM_CACHE_WARNINGS)

# Print a warning when an improperly aligned redistribution is performed, 
# i.e., if an unnecessary permutation communication stage must take place
option(ELEM_UNALIGNED_WARNINGS 
       "Warn when performing unaligned redistributions" OFF)
mark_as_advanced(ELEM_UNALIGNED_WARNINGS)

# Print a warning if an opportunity was missed to implement a redistribution
# approach specifically for vectors (instead of matrices)
option(ELEM_VECTOR_WARNINGS 
       "Warn when vector redistribution chances are missed" OFF)
mark_as_advanced(ELEM_VECTOR_WARNINGS)

################################################################################
# Significant command-line variable definitions                                #
################################################################################

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

# TODO: Add discussion of "MATH_LIBS"

# TODO: Add discussion of "REFERENCE_ROOT"

# TODO: Check to see if there are any others...

################################################################################
# Build logic starts here                                                      #
################################################################################

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
  set(ELEM_HYBRID TRUE)
  set(ELEM_DEBUG TRUE)
elseif(CMAKE_BUILD_TYPE STREQUAL "HybridRelease")
  set(ELEM_HYBRID TRUE)
  set(ELEM_RELEASE TRUE)
elseif(CMAKE_BUILD_TYPE STREQUAL "PureDebug")
  set(ELEM_PURE TRUE)
  set(ELEM_DEBUG TRUE)
elseif(CMAKE_BUILD_TYPE STREQUAL "PureRelease")
  set(ELEM_PURE TRUE)
  set(ELEM_RELEASE TRUE)
else()
  message(FATAL_ERROR "CMAKE_BUILD_TYPE must be [Hybrid,Pure][Debug,Release]")
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} UPPER_BUILD_TYPE)

# Set the basic compile flags from the build type
if(NOT WIN32) # I forget why we need this guard
  set(LANGUAGES CXX C Fortran)
  foreach(LANG ${LANGUAGES})
    if(NOT ${LANG}_FLAGS)
      if(ELEM_DEBUG)
        # set(${LANG}_FLAGS "-O2 -g -Wunused-variable -Wunused-but-set-variable -Wunused-local-typedefs" CACHE STRING
        set(${LANG}_FLAGS "-O2 -g" CACHE STRING
          "${LANG} optimization/debug flags for ${UPPER_BUILD_TYPE} mode")
      else()
        set(${LANG}_FLAGS "-O3" CACHE STRING
          "${LANG} optimization flags for ${UPPER_BUILD_TYPE} mode")
      endif()
    endif()
  endforeach()
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(tests/Fortran)
include(tests/CXX)
include(tests/MPI)
include(tests/OpenMP)

include(tests/Qt5)
if(ELEM_HAVE_QT5)
  set(ELEM_HEADERS_PREMOC 
      "include/elemental/io/DisplayWindow-premoc.hpp;include/elemental/io/ComplexDisplayWindow-premoc.hpp")
  qt_wrap_cpp(elemental ELEM_MOC_SRC ${ELEM_HEADERS_PREMOC})
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
  set(ELEM_HAVE_VALGRIND TRUE)
endif()

include(tests/Math)

# Now append the extra flags
set(CXX_FLAGS "${CXX_FLAGS} ${EXTRA_FLAGS}")
set(CMAKE_CXX_FLAGS_${UPPER_BUILD_TYPE} ${CXX_FLAGS})

# Add the Parallel Multiple Relatively Robust Representations (PMRRR) project
add_subdirectory(external/pmrrr)

# Create the Elemental configuration header
configure_file(${PROJECT_SOURCE_DIR}/cmake/config.h.cmake
               ${PROJECT_BINARY_DIR}/include/elemental/config.h)
install(FILES ${PROJECT_BINARY_DIR}/include/elemental/config.h 
        DESTINATION include/elemental)

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
set(LAPACK_ADDONSTRING)
foreach(LIB ${LAPACK_ADDONS_LIBS})
  set(LAPACK_ADDONSTRING "${LAPACK_ADDONSTRING} ${LIB}")
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
set(QT5_LIBSTRING "-L${Qt5_LIBDIR} -lQt5Widgets -lQt5Gui -lQt5Core")
configure_file(${PROJECT_SOURCE_DIR}/cmake/ElemVars.cmake
               ${PROJECT_BINARY_DIR}/conf/ElemVars @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/conf/ElemVars DESTINATION conf)

# Grab all of the .c, .cpp, .h, and .hpp Elemental files
file(GLOB_RECURSE ELEM_CPP RELATIVE ${PROJECT_SOURCE_DIR} "src/*.c" "src/*.cpp")
file(GLOB_RECURSE ELEM_HEADERS RELATIVE ${PROJECT_SOURCE_DIR} 
  "include/*.h" "include/*.hpp")
set(ELEM_SRC "${ELEM_CPP};${ELEM_HEADERS};${ELEM_MOC_SRC}")

# The main library
set(LINK_LIBS ${MATH_LIBS} ${MPI_CXX_LIBRARIES})
if(ELEM_HAVE_QT5)
  set(LINK_LIBS ${LINK_LIBS} ${Qt5Widgets_LIBRARIES})  
endif()
add_library(elemental ${LIBRARY_TYPE} ${ELEM_SRC})
target_link_libraries(elemental pmrrr ${LINK_LIBS})
if(MPI_CXX_LINK_FLAGS AND MPI_C_LINK_FLAGS)
  set(MPI_LINK_FLAGS "${MPI_CXX_LINK_FLAGS} ${MPI_C_LINK_FLAGS}")
elseif(MPI_CXX_LINK_FLAGS)
  set(MPI_LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
elseif(MPI_C_LINK_FLAGS)
  set(MPI_LINK_FLAGS "${MPI_C_LINK_FLAGS}")
endif()
install(TARGETS elemental DESTINATION lib)

# Define the header-file preparation rules
set(PREPARED_HEADERS)
foreach(HEADER ${ELEM_HEADERS})
  add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/${HEADER}
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/${HEADER}
            ${PROJECT_BINARY_DIR}/${HEADER}
    DEPENDS "${PROJECT_SOURCE_DIR}/${HEADER}")
  list(APPEND PREPARED_HEADERS ${PROJECT_BINARY_DIR}/${HEADER})

  get_filename_component(HEADER_PATH ${HEADER} PATH)
  install(FILES ${PROJECT_BINARY_DIR}/${HEADER} DESTINATION ${HEADER_PATH})
endforeach()
add_custom_target(prepare_elemental_headers DEPENDS ${PREPARED_HEADERS})
add_dependencies(elemental prepare_elemental_headers)

# TODO: Ensure that source files depend on header files

# Make sure the Elemental headers can be found
include_directories("${PROJECT_BINARY_DIR}/include")

# Build the test drivers if necessary
if(ELEM_TESTS)
  set(TEST_DIR ${PROJECT_SOURCE_DIR}/tests)
  set(TEST_TYPES core blas-like lapack-like convex)
  foreach(TYPE ${TEST_TYPES})
    file(GLOB_RECURSE ${TYPE}_TESTS RELATIVE ${PROJECT_SOURCE_DIR}/tests/${TYPE}/ 
         "tests/${TYPE}/*.cpp")
    set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/tests/${TYPE}")
    foreach(TEST ${${TYPE}_TESTS})
      set(DRIVER ${TEST_DIR}/${TYPE}/${TEST})
      get_filename_component(TESTNAME ${TEST} NAME_WE)
      add_executable(tests-${TYPE}-${TESTNAME} ${DRIVER})
      set_source_files_properties(${DRIVER} PROPERTIES 
        OBJECT_DEPENDS "${PREPARED_HEADERS}")
      target_link_libraries(tests-${TYPE}-${TESTNAME} elemental)
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
if(ELEM_EXAMPLES)
  set(EXAMPLE_DIR ${PROJECT_SOURCE_DIR}/examples)
  set(EXAMPLE_TYPES convex core blas-like io lapack-like matrices)
  foreach(TYPE ${EXAMPLE_TYPES})
    file(GLOB_RECURSE ${TYPE}_EXAMPLES RELATIVE 
         ${PROJECT_SOURCE_DIR}/examples/${TYPE}/ "examples/${TYPE}/*.cpp")
    set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/examples/${TYPE}")
    foreach(EXAMPLE ${${TYPE}_EXAMPLES})
      set(DRIVER ${EXAMPLE_DIR}/${TYPE}/${EXAMPLE})
      get_filename_component(EXNAME ${EXAMPLE} NAME_WE)
      add_executable(examples-${TYPE}-${EXNAME} ${DRIVER})
      set_source_files_properties(${DRIVER} PROPERTIES 
        OBJECT_DEPENDS "${PREPARED_HEADERS}")
      target_link_libraries(examples-${TYPE}-${EXNAME} elemental)
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
if(ELEM_EXPERIMENTAL)
  set(EXPERIMENTAL_DIR ${PROJECT_SOURCE_DIR}/experimental)

  # Build the G3D example(s)
  set(G3D_EXPERS G3DGemm)
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/experimental/g3d")
  foreach(EXPER ${G3D_EXPERS})
    set(DRIVER ${EXPERIMENTAL_DIR}/g3d/${EXPER}.cpp)
    add_executable(experimental-g3d-${EXPER} ${DRIVER})
    set_source_files_properties(${DRIVER} PROPERTIES 
      OBJECT_DEPENDS "${PREPARED_HEADERS}")
    target_link_libraries(experimental-g3d-${EXPER} elemental)
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
  include(./cmake/ElemSub.cmake)
endif()
