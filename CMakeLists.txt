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
set(Elemental_VERSION_MINOR 59)

# Elemental supports four build modes:
#   HybridDebug, HybridRelease, PureDebug, and PureRelease
if(NOT (CMAKE_BUILD_TYPE MATCHES HybridDebug) AND
   NOT (CMAKE_BUILD_TYPE MATCHES HybridRelease) AND 
   NOT (CMAKE_BUILD_TYPE MATCHES PureDebug) AND
   NOT (CMAKE_BUILD_TYPE MATCHES PureRelease))
  message(FATAL_ERROR 
          "CMAKE_BUILD_TYPE must be one of: HybridDebug, HybridRelease, PureDebug, PureRelease")
endif(NOT (CMAKE_BUILD_TYPE MATCHES HybridDebug) AND
      NOT (CMAKE_BUILD_TYPE MATCHES HybridRelease) AND 
      NOT (CMAKE_BUILD_TYPE MATCHES PureDebug) AND
      NOT (CMAKE_BUILD_TYPE MATCHES PureRelease))
set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING 
    "Choose the build type, options are: HybridDebug, HybridRelease, PureDebug, and PureRelease" 
    FORCE)

option(BUILD_TESTS   "Build a collection of test executables"             ON)
option(BUILD_PMRRR   "Build and link against eigensolver PMRRR"           ON)

option(WITHOUT_COMPLEX "Avoid all usage of complex numbers" OFF)
if(NOT WITHOUT_COMPLEX)
  option(AVOID_COMPLEX_MPI "Avoid complex MPI routines for robustness" ON)
endif(NOT WITHOUT_COMPLEX)
option(CACHE_WARNINGS "Warns when using cache-unfriendly routines"          OFF)
option(UNALIGNED_WARNINGS "Warn when performing unaligned redistributions"  OFF)
option(VECTOR_WARNINGS "Warn when vector redistribution chances are missed" OFF)
option(ENABLE_ALL_DISTRIBUTED_DOT "Build every parallel dot product."       OFF)
option(POOL_MEMORY "Make Memory class accumulate memory until destruction"  OFF)
mark_as_advanced(AVOID_COMPLEX_MPI CACHE_WARNINGS UNALIGNED_WARNINGS 
                 VECTOR_WARNINGS ENABLE_ALL_DISTRIBUTED_DOT POOL_MEMORY)

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
  set(CXX_HYBRID_DEBUG_FLAGS "-O2 -g" CACHE STRING
      "CXX optimization/debug flags for hybrid OpenMP/MPI debug build")
endif(HYBRID_DEBUG)
if(HYBRID_RELEASE)
  set(CXX_HYBRID_RELEASE_FLAGS "-O3" CACHE STRING
      "CXX optimization flags for hybrid OpenMP/MPI release build")
endif(HYBRID_RELEASE)
if(PURE_DEBUG)
  set(CXX_PURE_DEBUG_FLAGS "-O2 -g" CACHE STRING
      "CXX optimization/debugging flags for pure MPI debug build")
endif(PURE_DEBUG)
if(PURE_RELEASE)
  set(CXX_PURE_RELEASE_FLAGS "-O3" CACHE STRING
      "CXX optimization flags for pure MPI release build")
endif(PURE_RELEASE)

if(NOT (CMAKE_CXX_COMPILER MATCHES "mpi"))
  message(WARNING "CXX compiler, ${CMAKE_CXX_COMPILER} does not appear to be an MPI compiler. Please manually set it via the environment variable CXX or through the '-DCMAKE_CXX_COMPILER' cmake flag. You should probably also clean this build folder first.")
endif(NOT (CMAKE_CXX_COMPILER MATCHES "mpi"))

if(HYBRID)
  if(OpenMP_CXX_FLAGS)
    set(HAVE_OPENMP TRUE)
    set(OpenMP_CXX_FLAGS ${OpenMP_CXX_FLAGS})
  else(OpenMP_CXX_FLAGS)
    find_package(OpenMP)
    if(OPENMP_FOUND)
      set(HAVE_OPENMP TRUE)
    else(OPENMP_FOUND)
      set(OpenMP_CXX_FLAGS "" CACHE STRING "OpenMP CXX FLAGS")
      message(FATAL_ERROR
              "${CMAKE_BUILD_TYPE} failed because OpenMP support was not found. If your system supports OpenMP, please manually specify OpenMP_CXX_FLAGS.")
    endif(OPENMP_FOUND)
  endif(OpenMP_CXX_FLAGS)
endif(HYBRID)

# Query the size of a void pointer in order to determine whether or not this is
# a 32-bit system
if(${CMAKE_SIZEOF_VOID_P} MATCHES 4)
  set(32_BIT_SYSTEM TRUE)
else(${CMAKE_SIZEOF_VOID_P} MATCHES 4)
  set(32_BIT_SYSTEM FALSE)    
endif(${CMAKE_SIZEOF_VOID_P} MATCHES 4)

# Since the CMake find_blas for some reason requires a Fortran compiler, we
# must build our own.
include(${PROJECT_SOURCE_DIR}/cmake/find_library_list.cmake)
if(PURE)
  if(PURE_MATH_LIBS)
    set(NEED_PURE_MATH FALSE)
  else(PURE_MATH_LIBS)
    set(NEED_PURE_MATH TRUE)
  endif(PURE_MATH_LIBS)
else(PURE)
  set(NEED_PURE_MATH FALSE)
endif(PURE)
if(HYBRID AND HAVE_OPENMP AND NOT HYBRID_MATH_LIBS)
  set(NEED_HYBRID_MATH TRUE)
else(HYBRID AND HAVE_OPENMP AND NOT HYBRID_MATH_LIBS)
  set(NEED_HYBRID_MATH FALSE)
endif(HYBRID AND HAVE_OPENMP AND NOT HYBRID_MATH_LIBS)
if(NEED_PURE_MATH OR NEED_HYBRID_MATH)
  set(PURE_MATH_DESC "Serial BLAS/LAPACK link flags")
  set(HYBRID_MATH_DESC "Threaded BLAS/LAPACK link flags")
  if(APPLE)
    if(NEED_PURE_MATH)
      set(PURE_MATH_LIBS "-framework vecLib" 
          CACHE STRING ${PURE_MATH_DESC})
    endif(NEED_PURE_MATH)
    if(NEED_HYBRID_MATH)
      set(HYBRID_MATH_LIBS "-framework vecLib" 
          CACHE STRING ${HYBRID_MATH_DESC})
    endif(NEED_HYBRID_MATH)
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(FOUND_MATH_LIBS FALSE)
    # Look for MKL first
    if(32_BIT_SYSTEM)
      set(MKL_SUBDIR lib/32)
    else(32_BIT_SYSTEM)
      set(MKL_SUBDIR lib/em64t)
    endif(32_BIT_SYSTEM)
    set(MKL_LIBS mkl_intel_lp64 mkl_intel_thread mkl_core guide iomp5)
    find_library_list(MKL_MATH_LIBS MKL_LIBS "" ${MKL_SUBDIR})
    if(MKL_MATH_LIBS)
      find_library(PTHREAD_LIB pthread DOC "Pthreads library")
      if(PTHREAD_LIB)
        if(NEED_PURE_MATH)
          set(PURE_MATH_LIBS ${MKL_MATH_LIBS} 
              CACHE STRING ${PURE_MATH_DESC})
          message(STATUS "Using MKL for serial BLAS/LAPACK.")
        endif(NEED_PURE_MATH)
        if(NEED_HYBRID_MATH)
          set(HYBRID_MATH_LIBS ${MKL_MATH_LIBS} 
              CACHE STRING ${HYBRID_MATH_DESC})
          message(STATUS "Using MKL for threaded BLAS/LAPACK.")
        endif(NEED_HYBRID_MATH)
        set(FOUND_MATH_LIBS TRUE)
      else(PTHREAD_LIB)
        message(STATUS "WARNING: Found MKL but could not find pthread lib")  
        if(NEED_PURE_MATH)
          set(PURE_MATH_LIBS "" 
              CACHE STRING ${PURE_MATH_DESC})
        endif(NEED_PURE_MATH)
        if(NEED_HYBRID_MATH)
          set(HYBRID_MATH_LIBS "" 
              CACHE STRING ${HYBRID_MATH_DESC})
        endif(NEED_HYBRID_MATH)
      endif(PTHREAD_LIB)
    endif(MKL_MATH_LIBS)
    if(NOT FOUND_MATH_LIBS)
      # Look for default BLAS and LAPACK
      set(REFERENCE_LIBS lapack blas)
      find_library_list(REF_MATH_LIBS REFERENCE_LIBS "" "")
      if(REF_MATH_LIBS)
        if(NEED_PURE_MATH)
          set(PURE_MATH_LIBS ${REF_MATH_LIBS}
              CACHE STRING ${PURE_MATH_DESC})
          message(STATUS "WARNING: Using reference serial BLAS/LAPACK.")
        endif(NEED_PURE_MATH)
        if(NEED_HYBRID_MATH)
          set(HYBRID_MATH_LIBS ${REF_MATH_LIBS}
              CACHE STRING ${HYBRID_MATH_DESC})
          message(STATUS "WARNING: Using reference threaded BLAS/LAPACK.")
        endif(NEED_HYBRID_MATH)
      else(REF_MATH_LIBS)
        if(NEED_PURE_MATH)
          set(PURE_MATH_LIBS "" 
              CACHE STRING ${PURE_MATH_DESC})
        endif(NEED_PURE_MATH)
        if(NEED_HYBRID_MATH)
          set(HYBRID_MATH_LIBS "" 
              CACHE STRING ${HYBRID_MATH_DESC})
        endif(NEED_HYBRID_MATH)
        message(FATAL_ERROR 
          "Could not find BLAS/LAPACK libs. Please manually specify PURE_MATH_LIBS and/or HYBRID_MATH_LIBS.")
      endif(REF_MATH_LIBS)
    endif(NOT FOUND_MATH_LIBS)
  endif(APPLE)
endif(NEED_PURE_MATH OR NEED_HYBRID_MATH)

# Attempt to detect the BLAS and LAPACK underscore conventions. 
# We currently only handle whether or not there is an underscore appended.
include(CheckFunctionExists)
if(PURE)
  set(CMAKE_REQUIRED_LIBRARIES ${PURE_MATH_LIBS})
  check_function_exists(daxpy HAVE_DAXPY)
  if(HAVE_DAXPY)
    set(PURE_BLAS_POST FALSE)
    set(PURE_BLAS_DEFS "")
  else(HAVE_DAXPY)
    check_function_exists(daxpy_ HAVE_DAXPY_POST)
    if(HAVE_DAXPY_POST)
      set(PURE_BLAS_POST TRUE)
      set(PURE_BLAS_DEFS "-DBLAS_POST")
    else(HAVE_DAXPY_POST)
      message(FATAL_ERROR "Could not determine BLAS format.")
    endif(HAVE_DAXPY_POST)
  endif(HAVE_DAXPY)
  check_function_exists(dpotrf HAVE_DPOTRF)
  if(HAVE_DPOTRF)
    set(PURE_LAPACK_POST FALSE)
    set(PURE_LAPACK_DEFS "")
  else(HAVE_DPOTRF)
    check_function_exists(dpotrf_ HAVE_DPOTRF_POST)
    if(HAVE_DPOTRF_POST)
      set(PURE_LAPACK_POST TRUE)
      set(PURE_LAPACK_DEFS "-DLAPACK_POST")
    else(HAVE_DPOTRF_POST)
      message(FATAL_ERROR "Could not determine LAPACK format.")
    endif(HAVE_DPOTRF_POST)
  endif(HAVE_DPOTRF)
endif(PURE)
if(HYBRID AND HAVE_OPENMP)
  set(CMAKE_REQUIRED_LIBRARIES ${HYBRID_MATH_LIBS})
  check_function_exists(daxpy HAVE_DAXPY)
  if(HAVE_DAXPY)
    set(HYBRID_BLAS_POST FALSE)
    set(HYBRID_BLAS_DEFS "")
  else(HAVE_DAXPY)
    check_function_exists(daxpy_ HAVE_DAXPY_POST)
    if(HAVE_DAXPY_POST)
      set(HYBRID_BLAS_POST TRUE)
      set(HYBRID_BLAS_DEFS "-DBLAS_POST")
    else(HAVE_DAXPY_POST)
      message(FATAL_ERROR "Could not determine BLAS format.")
    endif(HAVE_DAXPY_POST)
  endif(HAVE_DAXPY)
  check_function_exists(dpotrf HAVE_DPOTRF)
  if(HAVE_DPOTRF)
    set(HYBRID_LAPACK_POST FALSE)
    set(HYBRID_LAPACK_DEFS "")
  else(HAVE_DPOTRF)
    check_function_exists(dpotrf_ HAVE_DPOTRF_POST)
    if(HAVE_DPOTRF_POST)
      set(HYBRID_LAPACK_POST TRUE)
      set(HYBRID_LAPACK_DEFS "-DLAPACK_POST")
    else(HAVE_DPOTRF_POST)
      message(FATAL_ERROR "Could not determine LAPACK format.")
    endif(HAVE_DPOTRF_POST)
  endif(HAVE_DPOTRF)
endif(HYBRID AND HAVE_OPENMP)
if(PURE)
  set(BLAS_POST ${PURE_BLAS_POST})
  set(LAPACK_POST ${PURE_LAPACK_POST})
else(PURE)
  set(BLAS_POST ${HYBRID_BLAS_POST})
  set(LAPACK_POST ${HYBRID_LAPACK_POST})
endif(PURE)

# Attempt to detect whether or not we have MKL by looking for mkl_set_num_threads
set(CMAKE_REQUIRED_LIBRARIES ${PURE_MATH_LIBS})
check_function_exists(mkl_set_num_threads PURE_MATH_HAS_MKL)
set(CMAKE_REQUIRED_LIBRARIES ${HYBRID_MATH_LIBS})
check_function_exists(mkl_set_num_threads HYBRID_MATH_HAS_MKL)

if(BUILD_PMRRR)
  add_subdirectory(contrib/pmrrr)
endif(BUILD_PMRRR)

if(BUILD_PMRRR AND NOT FAILED_PMRRR)
  set(WITHOUT_PMRRR FALSE)
else(BUILD_PMRRR AND NOT FAILED_PMRRR)
  set(WITHOUT_PMRRR TRUE)
endif(BUILD_PMRRR AND NOT FAILED_PMRRR)

configure_file( ${CMAKE_SOURCE_DIR}/include/elemental/config.h.cmake
                ${CMAKE_BINARY_DIR}/include/elemental/config.h ) 

# Grab all of the .c, .cpp, .h, and .hpp Elemental files
file(GLOB_RECURSE ELEMENTAL_CPP RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
     "src/*.c" "src/*.cpp")
file(GLOB_RECURSE ELEMENTAL_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
     "include/*.h" "include/*.hpp")
if(WITHOUT_PMRRR)
  message(STATUS "Removing eigensolver source from Elemental")
  list(REMOVE_ITEM ELEMENTAL_CPP
       "src/lapack/HermitianEig/HermitianEig.cpp")
  list(REMOVE_ITEM ELEMENTAL_CPP
       "src/lapack/SkewHermitianEig/SkewHermitianEig.cpp")
  list(REMOVE_ITEM ELEMENTAL_CPP
       "src/lapack/GeneralizedHermitianEig/GeneralizedHermitianEig.cpp")
endif(WITHOUT_PMRRR)
set(ELEMENTAL_SRC "${ELEMENTAL_CPP};${ELEMENTAL_HEADERS}")

# Create the lists of test executables and configure the compile flags
# for the various build modes
if(BUILD_TESTS)
  set(TEST_DIR ${PROJECT_SOURCE_DIR}/test)

  set(CORE_TESTS DifferentGrids DistMatrix Matrix)
  set(BLAS_TESTS Gemm Hemm Her2k Herk Symm Symv Syr2k Syrk Trmm Trsm)
  set(LAPACK_TESTS Chol Hegst LU QR Tridiag Trinv UT)

  set(C_TESTS BasicTest)

  set(CMAKE_CXX_FLAGS_HYBRIDDEBUG
      "${CXX_HYBRID_DEBUG_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS_HYBRIDRELEASE
      "${CXX_HYBRID_RELEASE_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS_PUREDEBUG 
      "${CXX_PURE_DEBUG_FLAGS}")
  set(CMAKE_CXX_FLAGS_PURERELEASE
      "${CXX_PURE_RELEASE_FLAGS}")
  if(BUILD_PMRRR AND NOT FAILED_PMRRR)
    list(APPEND LAPACK_TESTS GeneralizedHermitianEig)
    list(APPEND LAPACK_TESTS HermitianEig)
  endif(BUILD_PMRRR AND NOT FAILED_PMRRR)
endif(BUILD_TESTS)

# Set the math libraries based upon the build configuration
if(PURE)
  set(MATH_LIBS ${PURE_MATH_LIBS})
else(PURE)
  set(MATH_LIBS ${HYBRID_MATH_LIBS})
endif(PURE)

# Create a dummy library in order to be able to force the math libraries
# to be linked last
add_library(cmake-dummy-lib STATIC cmake/CMakeDummyFunction.cpp)
target_link_libraries(cmake-dummy-lib ${MATH_LIBS})

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
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/cxx/core")
  foreach(TEST ${CORE_TESTS})
    add_executable(cxx-core-${TEST} ${TEST_DIR}/cxx/core/${TEST}.cpp)
    target_link_libraries(cxx-core-${TEST} elemental cmake-dummy-lib)
    set_target_properties(cxx-core-${TEST} 
                          PROPERTIES OUTPUT_NAME ${TEST}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
  endforeach(TEST)
  # Build the blas tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/cxx/blas")
  foreach(TEST ${BLAS_TESTS})
    add_executable(cxx-blas-${TEST} ${TEST_DIR}/cxx/blas/${TEST}.cpp)
    target_link_libraries(cxx-blas-${TEST} elemental cmake-dummy-lib)
    set_target_properties(cxx-blas-${TEST} 
                          PROPERTIES OUTPUT_NAME ${TEST}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
  endforeach(TEST)
  # Build the lapack tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/cxx/lapack")
  foreach(TEST ${LAPACK_TESTS})
    add_executable(cxx-lapack-${TEST} ${TEST_DIR}/cxx/lapack/${TEST}.cpp)
    target_link_libraries(cxx-lapack-${TEST} elemental cmake-dummy-lib)
    set_target_properties(cxx-lapack-${TEST} 
                          PROPERTIES OUTPUT_NAME ${TEST}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
  endforeach(TEST)
  # Build the C tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/c")
  foreach(TEST ${C_TESTS})
    add_executable(c-${TEST} ${TEST_DIR}/c/${TEST}.c)
    target_link_libraries(c-${TEST} elemental cmake-dummy-lib)
    set_target_properties(c-${TEST}
                          PROPERTIES OUTPUT_NAME ${TEST}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
  endforeach(TEST)
endif(BUILD_TESTS)

