#
#  Copyright (c) 2009-2011, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#
#   - Redistributions of source code must retain the above copyright notice,
#     this list of conditions and the following disclaimer.
#
#   - Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#
#   - Neither the name of the owner nor the names of its contributors
#     may be used to endorse or promote products derived from this software
#     without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#  POSSIBILITY OF SUCH DAMAGE.
#
cmake_minimum_required(VERSION 2.8)
project(Elemental)

set(Elemental_VERSION_MAJOR 0)
set(Elemental_VERSION_MINOR 63)

# Elemental supports four build modes:
#   HybridDebug, HybridRelease, PureDebug, and PureRelease
if(NOT (CMAKE_BUILD_TYPE MATCHES HybridDebug) AND
   NOT (CMAKE_BUILD_TYPE MATCHES HybridRelease) AND 
   NOT (CMAKE_BUILD_TYPE MATCHES PureDebug) AND
   NOT (CMAKE_BUILD_TYPE MATCHES PureRelease))
  message(
    FATAL_ERROR 
    "CMAKE_BUILD_TYPE must be one of: HybridDebug, HybridRelease, PureDebug, PureRelease"
  )
endif(NOT (CMAKE_BUILD_TYPE MATCHES HybridDebug) AND
      NOT (CMAKE_BUILD_TYPE MATCHES HybridRelease) AND 
      NOT (CMAKE_BUILD_TYPE MATCHES PureDebug) AND
      NOT (CMAKE_BUILD_TYPE MATCHES PureRelease))
set(
  CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING 
  "Choose the build type, options are: HybridDebug, HybridRelease, PureDebug, and PureRelease" 
  FORCE
)

option(BUILD_TESTS   "Build a collection of test executables"   ON)
option(BUILD_EXAMPLES "Build the example(s)"                    ON)
option(BUILD_PMRRR   "Build and link against eigensolver PMRRR" ON)

option(WITHOUT_COMPLEX "Avoid all usage of complex numbers" OFF)
if(NOT WITHOUT_COMPLEX)
  option(AVOID_COMPLEX_MPI "Avoid complex MPI routines for robustness" ON)
endif(NOT WITHOUT_COMPLEX)
option(CACHE_WARNINGS "Warns when using cache-unfriendly routines"          OFF)
option(UNALIGNED_WARNINGS "Warn when performing unaligned redistributions"  OFF)
option(VECTOR_WARNINGS "Warn when vector redistribution chances are missed" OFF)
option(ENABLE_ALL_DISTRIBUTED_DOT "Build every parallel dot product."       OFF)
option(POOL_MEMORY "Make Memory class accumulate memory until destruction"  OFF)
mark_as_advanced(
  AVOID_COMPLEX_MPI CACHE_WARNINGS UNALIGNED_WARNINGS VECTOR_WARNINGS 
  ENABLE_ALL_DISTRIBUTED_DOT POOL_MEMORY
)

# Bug fix for IBM's compiler miscompiling OpenMP a[i] += alpha b[i] statements
# for complex datatypes
option(AVOID_OMP_FMA "Avoid a bug in the IBM compilers." ON)
mark_as_advanced(AVOID_OMP_FMA)

# Choose whether or not to build the prototype for the C interface
option(BUILD_C_INTERFACE "Build the prototype for the C interface" OFF)
mark_as_advanced(BUILD_C_INTERFACE)

if(CMAKE_BUILD_TYPE STREQUAL "HybridDebug")
  set(HYBRID TRUE)
  set(DEBUG TRUE)
  set(HYBRID_DEBUG TRUE)
endif(CMAKE_BUILD_TYPE STREQUAL "HybridDebug")
if(CMAKE_BUILD_TYPE STREQUAL "HybridRelease")
  set(HYBRID TRUE)
  set(RELEASE TRUE)
  set(HYBRID_RELEASE TRUE)
endif(CMAKE_BUILD_TYPE STREQUAL "HybridRelease")
if(CMAKE_BUILD_TYPE STREQUAL "PureDebug")
  set(PURE TRUE)
  set(DEBUG TRUE)
  set(PURE_DEBUG TRUE)
endif(CMAKE_BUILD_TYPE STREQUAL "PureDebug")
if(CMAKE_BUILD_TYPE STREQUAL "PureRelease")
  set(PURE TRUE)
  set(RELEASE TRUE)
  set(PURE_RELEASE TRUE)   
endif(CMAKE_BUILD_TYPE STREQUAL "PureRelease")

if(HYBRID_DEBUG)
  set(
    CXX_HYBRID_DEBUG_FLAGS "-O2 -g" CACHE STRING
    "CXX optimization/debug flags for hybrid OpenMP/MPI debug build"
  )
endif(HYBRID_DEBUG)
if(HYBRID_RELEASE)
  set(
    CXX_HYBRID_RELEASE_FLAGS "-O3" CACHE STRING
    "CXX optimization flags for hybrid OpenMP/MPI release build"
  )
endif(HYBRID_RELEASE)
if(PURE_DEBUG)
  set(
    CXX_PURE_DEBUG_FLAGS "-O2 -g" CACHE STRING
    "CXX optimization/debugging flags for pure MPI debug build"
  )
endif(PURE_DEBUG)
if(PURE_RELEASE)
  set(
    CXX_PURE_RELEASE_FLAGS "-O3" CACHE STRING
    "CXX optimization flags for pure MPI release build"
  )
endif(PURE_RELEASE)

if(NOT MPI_LINK_LIBS)
  if(NOT (CMAKE_CXX_COMPILER MATCHES "mpi"))
    message(
      WARNING 
      "CXX compiler, ${CMAKE_CXX_COMPILER} does not appear to be an MPI compiler. Please manually set it via the environment variable CXX or through the '-DCMAKE_CXX_COMPILER' cmake flag. You should probably also clean this build folder first."
    )
  endif(NOT (CMAKE_CXX_COMPILER MATCHES "mpi"))
endif(NOT MPI_LINK_LIBS)

if(HYBRID)
  if(OpenMP_CXX_FLAGS)
    set(HAVE_OPENMP TRUE)
    set(OpenMP_CXX_FLAGS ${OpenMP_CXX_FLAGS})
    message(STATUS "Using user-defined OpenMP_CXX_FLAGS=${OpenMP_CXX_FLAGS}")
  else(OpenMP_CXX_FLAGS)
    find_package(OpenMP)
    if(OPENMP_FOUND)
      set(HAVE_OPENMP TRUE)
    else(OPENMP_FOUND)
      set(OpenMP_CXX_FLAGS "" CACHE STRING "OpenMP CXX FLAGS")
      message(
        FATAL_ERROR
        "${CMAKE_BUILD_TYPE} failed because OpenMP support was not found. If your system supports OpenMP, please manually specify OpenMP_CXX_FLAGS."
      )
    endif(OPENMP_FOUND)
  endif(OpenMP_CXX_FLAGS)
endif(HYBRID)

set(
  CMAKE_CXX_FLAGS_HYBRIDDEBUG 
  "${CXX_HYBRID_DEBUG_FLAGS} ${OpenMP_CXX_FLAGS}"
)
set(
  CMAKE_CXX_FLAGS_HYBRIDRELEASE 
  "${CXX_HYBRID_RELEASE_FLAGS} ${OpenMP_CXX_FLAGS}"
)
set(CMAKE_CXX_FLAGS_PUREDEBUG "${CXX_PURE_DEBUG_FLAGS}")
set(CMAKE_CXX_FLAGS_PURERELEASE "${CXX_PURE_RELEASE_FLAGS}")

# Query the size of a void pointer in order to determine whether or not this is
# a 32-bit system
if(${CMAKE_SIZEOF_VOID_P} MATCHES 4)
  set(32_BIT_SYSTEM TRUE)
else(${CMAKE_SIZEOF_VOID_P} MATCHES 4)
  set(32_BIT_SYSTEM FALSE)    
endif(${CMAKE_SIZEOF_VOID_P} MATCHES 4)

# Since the CMake find_blas for some reason requires a Fortran compiler, we
# must build our own.
if(MATH_LIBS)
  set(NEED_MATH FALSE)
  message(STATUS "Using user-defined MATH_LIBS=${MATH_LIBS}")
else(MATH_LIBS)
  set(NEED_MATH TRUE)
endif(MATH_LIBS)
if(NEED_MATH)
  if(PURE)
    set(MATH_DESC "Unthreaded BLAS/LAPACK link flags")
  else(PURE)
    set(MATH_DESC "Threaded BLAS/LAPACK link flags")
  endif(PURE)
  if(APPLE)
    if(NEED_MATH)
      set(MATH_LIBS "-framework vecLib" CACHE STRING ${MATH_DESC})
      message(STATUS "Using Apple vecLib framework.")
    endif(NEED_MATH)
  else(APPLE)
    set(FOUND_MATH_LIBS FALSE)
    # Look for MKL first
    if(32_BIT_SYSTEM)
      set(INTEL_ARCH_SUBDIRS lib/32)
      set(MKL_ARCH_SUBDIRS lib/32)
    else(32_BIT_SYSTEM)
      set(INTEL_ARCH_SUBDIRS lib/em64t lib/intel64)
      set(MKL_ARCH_SUBDIRS lib/em64t lib/intel64 mkl/lib/em64t)
    endif(32_BIT_SYSTEM)
    if(INTEL_ROOT)
      message(
        STATUS 
        "Will also search user-defined INTEL_ROOT=${INTEL_ROOT} for MKL"
      )
    endif(INTEL_ROOT)
    if(MKL_ROOT)
      message(
        STATUS
        "Will also search user-defined MKL_ROOT=${MKL_ROOT} for MKL"
      )
    endif(MKL_ROOT)
    # There is a cyclic dependency between MKL_CORE and 
    # MKL_SEQUENTIAL/MKL_INTEL_THREAD with some versions of MKL
    find_library(CMATH_LIB m DOC "C math library")
    find_library(PTHREAD_LIB pthread DOC "Pthreads library")
    find_library(
      IOMP5_LIB
      NAMES libiomp5.a iomp5.lib 
      PATHS ${INTEL_ROOT} ${MKL_ROOT}
      PATH_SUFFIXES ${INTEL_ARCH_SUBDIRS}
    )
    find_library(
      GUIDE_LIB
      NAMES libguide.a guide.lib 
      PATHS ${INTEL_ROOT} ${MKL_ROOT}
      PATH_SUFFIXES ${INTEL_ARCH_SUBDIRS}
    )
    find_library(
      MKL_CORE_LIB
      NAMES libmkl_core.a mkl_core.lib 
      PATHS ${INTEL_ROOT} ${MKL_ROOT}
      PATH_SUFFIXES ${MKL_ARCH_SUBDIRS}
    )
    find_library(
      MKL_INTEL_LIB
      NAMES libmkl_intel_lp64.a mkl_intel_lp64.lib 
      PATHS ${INTEL_ROOT} ${MKL_ROOT}
      PATH_SUFFIXES ${MKL_ARCH_SUBDIRS}
    )
    if(PURE)
      find_library(
        MKL_SEQUENTIAL_LIB
        NAMES libmkl_sequential.a mkl_sequential.lib
        PATHS ${INTEL_ROOT} ${MKL_ROOT}
        PATH_SUFFIXES ${MKL_ARCH_SUBDIRS}
      )
      set(
        MKL_REQUIRED MKL_INTEL MKL_SEQUENTIAL MKL_CORE MKL_SEQUENTIAL 
        GUIDE IOMP5 PTHREAD CMATH
      ) 
    else(PURE)
      find_library(
        MKL_INTEL_THREAD_LIB
        NAMES libmkl_intel_thread.a mkl_intel_thread.lib
        PATH_SUFFIXES ${MKL_ARCH_SUBDIRS}
      )
      set(
        MKL_REQUIRED MKL_INTEL MKL_INTEL_THREAD MKL_CORE 
        MKL_INTEL_THREAD GUIDE IOMP5 PTHREAD CMATH
      )
    endif(PURE)
    set(MKL_FOUND TRUE)
    set(MATH_LIBS "")
    foreach(NAME ${MKL_REQUIRED})
      if(${NAME}_LIB)
        message(STATUS "Found ${NAME}_LIB: ${${NAME}_LIB}")
        list(APPEND MATH_LIBS ${${NAME}_LIB})
      else(${NAME}_LIB)
        message(STATUS "Could not find ${NAME}_LIB") 
        set(MATH_LIBS "")
        set(MKL_FOUND FALSE)
      endif(${NAME}_LIB)
    endforeach(NAME)
    if(MKL_FOUND)
      message(STATUS "Using MKL's math libraries.")
      set(NEED_MATH FALSE)
    endif(MKL_FOUND)
    if(NEED_MATH)
      # Look for default BLAS and LAPACK
      if(REFERENCE_ROOT)
        message(
          STATUS 
          "Will also search REFERENCE_ROOT=${REFERENCE_ROOT} for reference math libs"
        )
      endif(REFERENCE_ROOT)
      set(REFERENCE_REQUIRED LAPACK BLAS)
      find_library(BLAS_LIB NAMES blas PATHS ${REFERENCE_ROOT})
      find_library(LAPACK_LIB NAMES lapack reflapack PATHS ${REFERENCE_ROOT}) 
      set(REFERENCE_FOUND TRUE)
      foreach(NAME ${REFERENCE_REQUIRED})
        if(${NAME}_LIB)
          message(STATUS "Found ${NAME}_LIB: ${${NAME}_LIB}")
          list(APPEND MATH_LIBS ${${NAME}_LIB})
        else(${NAME}_LIB)
          message(STATUS "Could not find ${NAME}_LIB")
          set(MATH_LIBS "")
          set(REFERENCE_FOUND FALSE)
        endif(${NAME}_LIB)
      endforeach(NAME)
      if(REFERENCE_FOUND)
        message(STATUS "WARNING: Using reference BLAS/LAPACK.")
      else(REFERENCE_FOUND) 
        message(
          FATAL_ERROR 
          "Could not find BLAS/LAPACK libraries. Please manually specify MATH_LIBS"
        )
      endif(REFERENCE_FOUND)
    endif(NEED_MATH)
  endif(APPLE)
endif(NEED_MATH)

# Attempt to detect the BLAS and LAPACK underscore conventions. 
# We currently only handle whether or not there is an underscore appended.
include(CheckFunctionExists)
set(CMAKE_REQUIRED_LIBRARIES ${MATH_LIBS})
check_function_exists(daxpy HAVE_DAXPY)
if(HAVE_DAXPY)
  set(BLAS_POST FALSE)
  set(BLAS_DEFS "")
else(HAVE_DAXPY)
  check_function_exists(daxpy_ HAVE_DAXPY_POST)
  if(HAVE_DAXPY_POST)
    set(BLAS_POST TRUE)
    set(BLAS_DEFS "-DBLAS_POST")
  else(HAVE_DAXPY_POST)
    message(FATAL_ERROR "Could not determine BLAS format.")
  endif(HAVE_DAXPY_POST)
endif(HAVE_DAXPY)
check_function_exists(dpotrf HAVE_DPOTRF)
if(HAVE_DPOTRF)
  set(LAPACK_POST FALSE)
  set(LAPACK_DEFS "")
else(HAVE_DPOTRF)
  check_function_exists(dpotrf_ HAVE_DPOTRF_POST)
  if(HAVE_DPOTRF_POST)
    set(LAPACK_POST TRUE)
    set(LAPACK_DEFS "-DLAPACK_POST")
  else(HAVE_DPOTRF_POST)
    message(FATAL_ERROR "Could not determine LAPACK format.")
  endif(HAVE_DPOTRF_POST)
endif(HAVE_DPOTRF)

# Attempt to detect whether or not we have MKL by looking for 
# mkl_set_num_threads
set(CMAKE_REQUIRED_LIBRARIES ${MATH_LIBS})
check_function_exists(mkl_set_num_threads MATH_HAS_MKL)

if(BUILD_PMRRR)
  add_subdirectory(contrib/pmrrr)
endif(BUILD_PMRRR)

if(BUILD_PMRRR AND NOT FAILED_PMRRR)
  set(WITHOUT_PMRRR FALSE)
else(BUILD_PMRRR AND NOT FAILED_PMRRR)
  set(WITHOUT_PMRRR TRUE)
endif(BUILD_PMRRR AND NOT FAILED_PMRRR)

configure_file( 
  ${CMAKE_SOURCE_DIR}/include/elemental/config.h.cmake
  ${CMAKE_BINARY_DIR}/include/elemental/config.h 
) 

# Grab all of the .c, .cpp, .h, and .hpp Elemental files
file(
  GLOB_RECURSE ELEMENTAL_CPP RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
  "src/*.c" "src/*.cpp"
)
file(
  GLOB_RECURSE ELEMENTAL_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
  "include/*.h" "include/*.hpp"
)
if(WITHOUT_PMRRR)
  message(STATUS "Removing eigensolver source from Elemental")
  list(
    REMOVE_ITEM ELEMENTAL_CPP
    "src/lapack/HermitianEig/HermitianEig.cpp"
  )
  list(
    REMOVE_ITEM ELEMENTAL_CPP
   "src/lapack/SkewHermitianEig/SkewHermitianEig.cpp"
  )
  list(
    REMOVE_ITEM ELEMENTAL_CPP
    "src/lapack/GeneralizedHermitianEig/GeneralizedHermitianEig.cpp"
  )
endif(WITHOUT_PMRRR)
if(NOT BUILD_C_INTERFACE)
    message(STATUS "Removing prototype C interface from Elemental")
    list(REMOVE_ITEM ELEMENTAL_CPP "src/exports/c_interface.cpp")
    list(REMOVE_ITEM ELEMENTAL_HEADERS "include/elemental/exports/c.h")
endif(NOT BUILD_C_INTERFACE)
set(ELEMENTAL_SRC "${ELEMENTAL_CPP};${ELEMENTAL_HEADERS}")

# Create the lists of test executables 
if(BUILD_TESTS)
  set(TEST_DIR ${PROJECT_SOURCE_DIR}/test)

  set(CORE_TESTS AxpyInterface DifferentGrids DistMatrix Matrix)
  set(BASIC_TESTS Gemm Hemm Her2k Herk Symm Symv Syr2k Syrk Trmm Trsm)
  set(ADVANCED_TESTS Chol Hegst LU QR Tridiag Trinv UT)

  if(BUILD_C_INTERFACE)
    set(C_TESTS BasicTest)
  endif(BUILD_C_INTERFACE)

  if(BUILD_PMRRR AND NOT FAILED_PMRRR)
    list(APPEND ADVANCED_TESTS GeneralizedHermitianEig)
    list(APPEND ADVANCED_TESTS HermitianEig)
  endif(BUILD_PMRRR AND NOT FAILED_PMRRR)
endif(BUILD_TESTS)

# Create the lists of example executables. There is currently only one.
if(BUILD_EXAMPLES)
  set(EXAMPLE_DIR ${PROJECT_SOURCE_DIR}/examples)

  set(ADVANCED_EXAMPLES)

  if(BUILD_PMRRR AND NOT FAILED_PMRRR)
    if(NOT WITHOUT_COMPLEX)
      list(APPEND ADVANCED_EXAMPLES HermitianEig)
    endif(NOT WITHOUT_COMPLEX)
    list(APPEND ADVANCED_EXAMPLES SymmetricEig)
  endif(BUILD_PMRRR AND NOT FAILED_PMRRR)
endif(BUILD_EXAMPLES)

# Create a dummy library in order to be able to force the math (and potentially
# MPI libraries) to be linked last
add_library(cmake-dummy-lib STATIC cmake/CMakeDummyFunction.cpp)
target_link_libraries(cmake-dummy-lib ${MATH_LIBS} ${MPI_LINK_LIBS})

# If we specified an MPI include directory, make sure it's in our include path
if(MPI_INCLUDE_DIR)
  include_directories("${MPI_INCLUDE_DIR}")
endif(MPI_INCLUDE_DIR)

# Copy the headers into the build directory
foreach(HEADER ${ELEMENTAL_HEADERS})
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${HEADER}
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${HEADER}
            ${CMAKE_CURRENT_BINARY_DIR}/${HEADER}
  )
endforeach(HEADER)

# Make sure the Elemental headers can be found
include_directories("${PROJECT_BINARY_DIR}/include")

# Build the main library
if(BUILD_PMRRR AND NOT FAILED_PMRRR)
  add_library(elemental STATIC ${ELEMENTAL_SRC})
  target_link_libraries(elemental pmrrr)
else(BUILD_PMRRR AND NOT FAILED_PMRRR)
  add_library(elemental STATIC ${ELEMENTAL_SRC})
endif(BUILD_PMRRR AND NOT FAILED_PMRRR)

# Build the test drivers if necessary
if(BUILD_TESTS)
  # Build the core tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/tests/cxx/core")
  foreach(TEST ${CORE_TESTS})
    add_executable(tests-cxx-core-${TEST} ${TEST_DIR}/cxx/core/${TEST}.cpp)
    target_link_libraries(tests-cxx-core-${TEST} elemental cmake-dummy-lib)
    set_target_properties(
      tests-cxx-core-${TEST} 
      PROPERTIES OUTPUT_NAME ${TEST}
      RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
    )
  endforeach(TEST)
  # Build the basic tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/tests/cxx/basic")
  foreach(TEST ${BASIC_TESTS})
    add_executable(tests-cxx-basic-${TEST} ${TEST_DIR}/cxx/basic/${TEST}.cpp)
    target_link_libraries(tests-cxx-basic-${TEST} elemental cmake-dummy-lib)
    set_target_properties(
      tests-cxx-basic-${TEST} 
      PROPERTIES OUTPUT_NAME ${TEST}
      RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
    )
  endforeach(TEST)
  # Build the advanced tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/tests/cxx/advanced")
  foreach(TEST ${ADVANCED_TESTS})
    add_executable(
      tests-cxx-advanced-${TEST} 
      ${TEST_DIR}/cxx/advanced/${TEST}.cpp
    )
    target_link_libraries(tests-cxx-advanced-${TEST} elemental cmake-dummy-lib)
    set_target_properties(
      tests-cxx-advanced-${TEST} 
      PROPERTIES OUTPUT_NAME ${TEST}
      RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
    )
  endforeach(TEST)
  # Build the C tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/tests/c")
  foreach(TEST ${C_TESTS})
    add_executable(tests-c-${TEST} ${TEST_DIR}/c/${TEST}.c)
    target_link_libraries(tests-c-${TEST} elemental cmake-dummy-lib)
    set_target_properties(
      tests-c-${TEST}
      PROPERTIES OUTPUT_NAME ${TEST}
      RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
    )
  endforeach(TEST)
endif(BUILD_TESTS)

# Build the example drivers if necessary
if(BUILD_EXAMPLES)
  # Build the advanced test(s)
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/examples/cxx/advanced")
  foreach(EXAMPLE ${ADVANCED_EXAMPLES})
    add_executable(
      examples-cxx-advanced-${EXAMPLE} 
      ${EXAMPLE_DIR}/cxx/advanced/${EXAMPLE}.cpp
    )
    target_link_libraries(
      examples-cxx-advanced-${EXAMPLE} 
      elemental cmake-dummy-lib
    )
    set_target_properties(
      examples-cxx-advanced-${EXAMPLE} 
      PROPERTIES OUTPUT_NAME ${EXAMPLE}
      RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
    )
  endforeach(EXAMPLE)
endif(BUILD_EXAMPLES)

