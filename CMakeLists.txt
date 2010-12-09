#
#  Copyright (c) 2009-2010, Jack Poulson
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
cmake_minimum_required(VERSION 2.6)
project(Elemental)
set(Elemental_VERSION_MAJOR 0)
set(Elemental_VERSION_MINOR 53)

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
option(TIMING "If defined, some routines will accumulate statistics"        ON)
option(CACHE_WARNINGS "Warns when using cache-unfriendly routines"          OFF)
option(UNALIGNED_WARNINGS "Warn when performing unaligned redistributions"  OFF)
option(VECTOR_WARNINGS "Warn when vector redistribution chances are missed" OFF)
option(ENABLE_ALL_DISTRIBUTED_DOT "Build every parallel dot product."       OFF)
option(POOL_MEMORY "Make Memory class accumulate memory until destruction"  OFF)
mark_as_advanced(AVOID_COMPLEX_MPI TIMING CACHE_WARNINGS UNALIGNED_WARNINGS 
                 VECTOR_WARNINGS ENABLE_ALL_DISTRIBUTED_DOT POOL_MEMORY)

if(CMAKE_BUILD_TYPE STREQUAL "HybridDebug")
  set(BUILD_HYBRID TRUE)
  set(BUILD_DEBUG TRUE)
  set(BUILD_HYBRID_DEBUG TRUE)
endif(CMAKE_BUILD_TYPE STREQUAL "HybridDebug")
if(CMAKE_BUILD_TYPE STREQUAL "HybridRelease")
  set(BUILD_HYBRID TRUE)
  set(BUILD_RELEASE TRUE)
  set(BUILD_HYBRID_RELEASE TRUE)
endif(CMAKE_BUILD_TYPE STREQUAL "HybridRelease")
if(CMAKE_BUILD_TYPE STREQUAL "PureDebug")
  set(BUILD_PURE TRUE)
  set(BUILD_DEBUG TRUE)
  set(BUILD_PURE_DEBUG TRUE)
endif(CMAKE_BUILD_TYPE STREQUAL "PureDebug")
if(CMAKE_BUILD_TYPE STREQUAL "PureRelease")
  set(BUILD_PURE TRUE)
  set(BUILD_RELEASE TRUE)
  set(BUILD_PURE_RELEASE TRUE)   
endif(CMAKE_BUILD_TYPE STREQUAL "PureRelease")

if(BUILD_HYBRID_DEBUG)
  set(CXX_HYBRID_DEBUG_FLAGS "-O2 -g" CACHE STRING
      "CXX optimization/debug flags for hybrid OpenMP/MPI debug build")
endif(BUILD_HYBRID_DEBUG)
if(BUILD_HYBRID_RELEASE)
  set(CXX_HYBRID_RELEASE_FLAGS "-O3" CACHE STRING
      "CXX optimization flags for hybrid OpenMP/MPI release build")
endif(BUILD_HYBRID_RELEASE)
if(BUILD_PURE_DEBUG)
  set(CXX_PURE_DEBUG_FLAGS "-O2 -g" CACHE STRING
      "CXX optimization/debugging flags for pure MPI debug build")
endif(BUILD_PURE_DEBUG)
if(BUILD_PURE_RELEASE)
  set(CXX_PURE_RELEASE_FLAGS "-O3" CACHE STRING
      "CXX optimization flags for pure MPI release build")
endif(BUILD_PURE_RELEASE)

if(NOT (CMAKE_CXX_COMPILER MATCHES "mpi"))
  message(WARNING "CXX compiler, ${CMAKE_CXX_COMPILER} does not appear to be an MPI compiler. Please manually set it via the environment variable CXX or through the '-DCMAKE_CXX_COMPILER' cmake flag. You should probably also clean this build folder first.")
endif(NOT (CMAKE_CXX_COMPILER MATCHES "mpi"))

if(BUILD_HYBRID)
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
endif(BUILD_HYBRID)

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
if(BUILD_PURE OR BUILD_PMRRR)
  if(PURE_BLAS_LAPACK_LIBS)
    set(NEED_PURE_BLAS_LAPACK FALSE)
  else(PURE_BLAS_LAPACK_LIBS)
    set(NEED_PURE_BLAS_LAPACK TRUE)
  endif(PURE_BLAS_LAPACK_LIBS)
else(BUILD_PURE OR BUILD_PMRRR)
  set(NEED_PURE_BLAS_LAPACK FALSE)
endif(BUILD_PURE OR BUILD_PMRRR)
if(BUILD_HYBRID AND HAVE_OPENMP AND NOT HYBRID_BLAS_LAPACK_LIBS)
  set(NEED_HYBRID_BLAS_LAPACK TRUE)
else(BUILD_HYBRID AND HAVE_OPENMP AND NOT HYBRID_BLAS_LAPACK_LIBS)
  set(NEED_HYBRID_BLAS_LAPACK FALSE)
endif(BUILD_HYBRID AND HAVE_OPENMP AND NOT HYBRID_BLAS_LAPACK_LIBS)
if(NEED_PURE_BLAS_LAPACK OR NEED_HYBRID_BLAS_LAPACK)
  set(PURE_BLAS_LAPACK_DESC "Serial BLAS/LAPACK link flags")
  set(HYBRID_BLAS_LAPACK_DESC "Threaded BLAS/LAPACK link flags")
  if(APPLE)
    if(NEED_PURE_BLAS_LAPACK)
      set(PURE_BLAS_LAPACK_LIBS "-framework vecLib" 
          CACHE STRING ${PURE_BLAS_LAPACK_DESC})
    endif(NEED_PURE_BLAS_LAPACK)
    if(NEED_HYBRID_BLAS_LAPACK)
      set(HYBRID_BLAS_LAPACK_LIBS "-framework vecLib" 
          CACHE STRING ${HYBRID_BLAS_LAPACK_DESC})
    endif(NEED_HYBRID_BLAS_LAPACK)
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(FOUND_BLAS_LAPACK_LIBS FALSE)
    # Look for MKL first
    if(32_BIT_SYSTEM)
      set(MKL_SUBDIR lib/32)
    else(32_BIT_SYSTEM)
      set(MKL_SUBDIR lib/em64t)
    endif(32_BIT_SYSTEM)
    set(MKL_LIBS mkl_intel_lp64 mkl_intel_thread mkl_core guide iomp5)
    find_library_list(MKL_BLAS_LAPACK_LIBS MKL_LIBS "" ${MKL_SUBDIR})
    if(MKL_BLAS_LAPACK_LIBS)
      find_library(PTHREAD_LIB pthread DOC "Pthreads library")
      if(PTHREAD_LIB)
        if(NEED_PURE_BLAS_LAPACK)
          set(PURE_BLAS_LAPACK_LIBS ${MKL_BLAS_LAPACK_LIBS} 
              CACHE STRING ${PURE_BLAS_LAPACK_DESC})
          message(STATUS "Using MKL for serial BLAS/LAPACK.")
        endif(NEED_PURE_BLAS_LAPACK)
        if(NEED_HYBRID_BLAS_LAPACK)
          set(HYBRID_BLAS_LAPACK_LIBS ${MKL_BLAS_LAPACK_LIBS} 
              CACHE STRING ${HYBRID_BLAS_LAPACK_DESC})
          message(STATUS "Using MKL for threaded BLAS/LAPACK.")
        endif(NEED_HYBRID_BLAS_LAPACK)
        set(FOUND_BLAS_LAPACK_LIBS TRUE)
      else(PTHREAD_LIB)
        message(STATUS "WARNING: Found MKL but could not find pthread lib")  
        if(NEED_PURE_BLAS_LAPACK)
          set(PURE_BLAS_LAPACK_LIBS "" 
              CACHE STRING ${PURE_BLAS_LAPACK_DESC})
        endif(NEED_PURE_BLAS_LAPACK)
        if(NEED_HYBRID_BLAS_LAPACK)
          set(HYBRID_BLAS_LAPACK_LIBS "" 
              CACHE STRING ${HYBRID_BLAS_LAPACK_DESC})
        endif(NEED_HYBRID_BLAS_LAPACK)
      endif(PTHREAD_LIB)
    endif(MKL_BLAS_LAPACK_LIBS)
    if(NOT FOUND_BLAS_LAPACK_LIBS)
      # Look for default BLAS and LAPACK
      set(REFERENCE_LIBS lapack blas)
      find_library_list(REF_BLAS_LAPACK_LIBS REFERENCE_LIBS "" "")
      if(REF_BLAS_LAPACK_LIBS)
        if(NEED_PURE_BLAS_LAPACK)
          set(PURE_BLAS_LAPACK_LIBS ${REF_BLAS_LAPACK_LIBS}
              CACHE STRING ${PURE_BLAS_LAPACK_DESC})
          message(STATUS "WARNING: Using reference serial BLAS/LAPACK.")
        endif(NEED_PURE_BLAS_LAPACK)
        if(NEED_HYBRID_BLAS_LAPACK)
          set(HYBRID_BLAS_LAPACK_LIBS ${REF_BLAS_LAPACK_LIBS}
              CACHE STRING ${HYBRID_BLAS_LAPACK_DESC})
          message(STATUS "WARNING: Using reference threaded BLAS/LAPACK.")
        endif(NEED_HYBRID_BLAS_LAPACK)
      else(REF_BLAS_LAPACK_LIBS)
        if(NEED_PURE_BLAS_LAPACK)
          set(PURE_BLAS_LAPACK_LIBS "" 
              CACHE STRING ${PURE_BLAS_LAPACK_DESC})
        endif(NEED_PURE_BLAS_LAPACK)
        if(NEED_HYBRID_BLAS_LAPACK)
          set(HYBRID_BLAS_LAPACK_LIBS "" 
              CACHE STRING ${HYBRID_BLAS_LAPACK_DESC})
        endif(NEED_HYBRID_BLAS_LAPACK)
        message(FATAL_ERROR 
          "Could not find BLAS/LAPACK libs. Please manually specify BLAS_LAPACK_LIBS.")
      endif(REF_BLAS_LAPACK_LIBS)
    endif(NOT FOUND_BLAS_LAPACK_LIBS)
  endif(APPLE)
endif(NEED_PURE_BLAS_LAPACK OR NEED_HYBRID_BLAS_LAPACK)

# Begin building CXX_BASE_FLAGS
set(CXX_BASE_FLAGS)
if(WITHOUT_COMPLEX)
  set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DWITHOUT_COMPLEX")
elseif(AVOID_COMPLEX_MPI)
  set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DAVOID_COMPLEX_MPI")
endif(WITHOUT_COMPLEX)
if(TIMING)
  set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DTIMING")
endif(TIMING)
if(CACHE_WARNINGS)
  set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DCACHE_WARNINGS")
endif(CACHE_WARNINGS)
if(UNALIGNED_WARNINGS)
  set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DUNALIGNED_WARNINGS")
endif(UNALIGNED_WARNINGS)
if(VECTOR_WARNINGS)
  set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DVECTOR_WARNINGS")
endif(VECTOR_WARNINGS)
if(ENABLE_ALL_DISTRIBUTED_DOT)
  set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DENABLE_ALL_DISTRIBUTED_DOT")
endif(ENABLE_ALL_DISTRIBUTED_DOT)
if(POOL_MEMORY)
  set(CXX_BASE_FLAGS "${CXX_BASE_FLAGS} -DPOOL_MEMORY")
endif(POOL_MEMORY)

# Attempt to detect the BLAS and LAPACK underscore conventions
include(CheckFunctionExists)
if(BUILD_PURE)
  set(CMAKE_REQUIRED_LIBRARIES ${PURE_BLAS_LAPACK_LIBS})
  check_function_exists(daxpy HAVE_DAXPY)
  if(HAVE_DAXPY)
    set(PURE_BLAS_DEFS "")
  else(HAVE_DAXPY)
    check_function_exists(_daxpy HAVE_DAXPY_PRE)
    if(HAVE_DAXPY_PRE)
      set(PURE_BLAS_DEFS "-DBLAS_PRE")
    else(HAVE_DAXPY_PRE)
      check_function_exists(daxpy_ HAVE_DAXPY_POST)
      if(HAVE_DAXPY_POST)
        set(PURE_BLAS_DEFS "-DBLAS_POST")
      else(HAVE_DAXPY_POST)
        check_function_exists(_daxpy_ HAVE_DAXPY_PRE_POST)
        if(HAVE_DAXPY_PRE_POST)
          set(PURE_BLAS_DEFS "-DBLAS_PRE -DBLAS_POST")
        else(HAVE_DAXPY_PRE_POST)
          message(FATAL_ERROR "Could not determine BLAS format.")
        endif(HAVE_DAXPY_PRE_POST)
      endif(HAVE_DAXPY_POST)
    endif(HAVE_DAXPY_PRE)
  endif(HAVE_DAXPY)
  check_function_exists(dpotrf HAVE_DPOTRF)
  if(HAVE_DPOTRF)
    set(PURE_LAPACK_DEFS "")
  else(HAVE_DPOTRF)
    check_function_exists(_dpotrf HAVE_DPOTRF_PRE)
    if(HAVE_DPOTRF_PRE)
      set(PURE_LAPACK_DEFS "-DLAPACK_PRE")
    else(HAVE_DPOTRF_PRE)
      check_function_exists(dpotrf_ HAVE_DPOTRF_POST)
      if(HAVE_DPOTRF_POST)
        set(PURE_LAPACK_DEFS "-DLAPACK_POST")
      else(HAVE_DPOTRF_POST)
        check_function_exists(_dpotrf_ HAVE_DPOTRF_PRE_POST)
        if(HAVE_DPOTRF_PRE_POST)
          set(PURE_LAPACK_DEFS "-DLAPACK_PRE -DLAPACK_POST")
        else(HAVE_DPOTRF_PRE_POST)
          message(FATAL_ERROR "Could not determine LAPACK format.")
        endif(HAVE_DPOTRF_PRE_POST)
      endif(HAVE_DPOTRF_POST)
    endif(HAVE_DPOTRF_PRE)
  endif(HAVE_DPOTRF)
endif(BUILD_PURE)
if(BUILD_HYBRID AND HAVE_OPENMP)
  set(CMAKE_REQUIRED_LIBRARIES ${HYBRID_BLAS_LAPACK_LIBS})
  check_function_exists(daxpy HAVE_DAXPY)
  if(HAVE_DAXPY)
    set(HYBRID_BLAS_DEFS "")
  else(HAVE_DAXPY)
    check_function_exists(_daxpy HAVE_DAXPY_PRE)
    if(HAVE_DAXPY_PRE)
      set(HYBRID_BLAS_DEFS "-DBLAS_PRE")
    else(HAVE_DAXPY_PRE)
      check_function_exists(daxpy_ HAVE_DAXPY_POST)
      if(HAVE_DAXPY_POST)
        set(HYBRID_BLAS_DEFS "-DBLAS_POST")
      else(HAVE_DAXPY_POST)
        check_function_exists(_daxpy_ HAVE_DAXPY_PRE_POST)
        if(HAVE_DAXPY_PRE_POST)
          set(HYBRID_BLAS_DEFS "-DBLAS_PRE -DBLAS_POST")
        else(HAVE_DAXPY_PRE_POST)
          message(FATAL_ERROR "Could not determine BLAS format.")
        endif(HAVE_DAXPY_PRE_POST)
      endif(HAVE_DAXPY_POST)
    endif(HAVE_DAXPY_PRE)
  endif(HAVE_DAXPY)
  check_function_exists(dpotrf HAVE_DPOTRF)
  if(HAVE_DPOTRF)
    set(HYBRID_LAPACK_DEFS "")
  else(HAVE_DPOTRF)
    check_function_exists(_dpotrf HAVE_DPOTRF_PRE)
    if(HAVE_DPOTRF_PRE)
      set(HYBRID_LAPACK_DEFS "-DLAPACK_PRE")
    else(HAVE_DPOTRF_PRE)
      check_function_exists(dpotrf_ HAVE_DPOTRF_POST)
      if(HAVE_DPOTRF_POST)
        set(HYBRID_LAPACK_DEFS "-DLAPACK_POST")
      else(HAVE_DPOTRF_POST)
        check_function_exists(_dpotrf_ HAVE_DPOTRF_PRE_POST)
        if(HAVE_DPOTRF_PRE_POST)
          set(HYBRID_LAPACK_DEFS "-DLAPACK_PRE -DLAPACK_POST")
        else(HAVE_DPOTRF_PRE_POST)
          message(FATAL_ERROR "Could not determine LAPACK format.")
        endif(HAVE_DPOTRF_PRE_POST)
      endif(HAVE_DPOTRF_POST)
    endif(HAVE_DPOTRF_PRE)
  endif(HAVE_DPOTRF)
endif(BUILD_HYBRID AND HAVE_OPENMP)

if(BUILD_PMRRR)
  add_subdirectory(contrib/pmrrr)
endif(BUILD_PMRRR)

# Grab all of the .cpp and .hpp Elemental files
file(GLOB_RECURSE ELEMENTAL_CPP RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
     "src/*.cpp")
file(GLOB_RECURSE ELEMENTAL_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
     "include/*.hpp")
set(ELEMENTAL_SRC "${ELEMENTAL_CPP};${ELEMENTAL_HEADERS}")
if(NOT BUILD_PMRRR OR FAILED_PMRRR)
  message(STATUS "Removing eigensolver source from Elemental")
  list(REMOVE_ITEM ELEMENTAL_SRC 
       "src/lapack/HermitianEig/HermitianEig.cpp")
  list(REMOVE_ITEM ELEMENTAL_SRC 
       "src/lapack/GeneralizedHermitianEig/GeneralizedHermitianEig.cpp")
endif(NOT BUILD_PMRRR OR FAILED_PMRRR)

# Create the lists of test executables and configure the compile flags
# for the various build modes
if(BUILD_TESTS)
  set(TEST_DIR ${PROJECT_SOURCE_DIR}/test)
  set(CORE_TESTS DifferentGrids DistMatrix Matrix)
  set(BLAS_TESTS Gemm Hemm Her2k Herk Symm Symv Syr2k Syrk Trmm Trsm)
  set(LAPACK_TESTS Chol Hegst LU QR Tridiag Trinv UT)

  set(CMAKE_CXX_FLAGS_HYBRIDDEBUG
      "${CXX_BASE_FLAGS} ${HYBRID_BLAS_DEFS} ${HYBRID_LAPACK_DEFS} ${CXX_HYBRID_DEBUG_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS_HYBRIDRELEASE
      "${CXX_BASE_FLAGS} ${HYBRID_BLAS_DEFS} ${HYBRID_LAPACK_DEFS} ${CXX_HYBRID_RELEASE_FLAGS} ${OpenMP_CXX_FLAGS} -DRELEASE")
  set(CMAKE_CXX_FLAGS_PUREDEBUG 
      "${CXX_BASE_FLAGS} ${PURE_BLAS_DEFS} ${PURE_LAPACK_DEFS} ${CXX_PURE_DEBUG_FLAGS}")
  set(CMAKE_CXX_FLAGS_PURERELEASE
      "${CXX_BASE_FLAGS} ${PURE_BLAS_DEFS} ${PURE_LAPACK_DEFS} ${CXX_PURE_RELEASE_FLAGS} -DRELEASE")
  if(BUILD_PMRRR AND NOT FAILED_PMRRR)
    list(APPEND LAPACK_TESTS GeneralizedHermitianEig)
    list(APPEND LAPACK_TESTS HermitianEig)
  else(BUILD_PMRRR AND NOT FAILED_PMRRR)
    set(CMAKE_CXX_FLAGS_HYBRIDDEBUG 
        "${CMAKE_CXX_FLAGS_HYBRIDDEBUG} -DWITHOUT_PMRRR")
    set(CMAKE_CXX_FLAGS_HYBRIDRELEASE 
        "${CMAKE_CXX_FLAGS_HYBRIDRELEASE} -DWITHOUT_PMRRR")
    set(CMAKE_CXX_FLAGS_PUREDEBUG 
        "${CMAKE_CXX_FLAGS_PUREDEBUG} -DWITHOUT_PMRRR")
    set(CMAKE_CXX_FLAGS_PURERELEASE 
        "${CMAKE_CXX_FLAGS_PURERELEASE} -DWITHOUT_PMRRR")
  endif(BUILD_PMRRR AND NOT FAILED_PMRRR)
endif(BUILD_TESTS)

# Make sure the Elemental headers can be found
include_directories("${PROJECT_SOURCE_DIR}/include")

# Unroll the BLAS/LAPACK libs into a string
if(BUILD_PURE)
  set(BLAS_LAPACK_FLAGS "")
  foreach(LIB ${PURE_BLAS_LAPACK_LIBS})
    set(BLAS_LAPACK_FLAGS "${BLAS_LAPACK_FLAGS} ${LIB}")
  endforeach(LIB)
else(BUILD_PURE)
  set(BLAS_LAPACK_FLAGS "")
  foreach(LIB ${HYBRID_BLAS_LAPACK_LIBS})
    set(BLAS_LAPACK_FLAGS "${BLAS_LAPACK_FLAGS} ${LIB}")
  endforeach(LIB)
endif(BUILD_PURE)

if(BUILD_PMRRR AND NOT FAILED_PMRRR)
  add_library(elemental STATIC ${ELEMENTAL_SRC})
  target_link_libraries(elemental pmrrr)
else(BUILD_PMRRR AND NOT FAILED_PMRRR)
  add_library(elemental STATIC ${ELEMENTAL_SRC})
endif(BUILD_PMRRR AND NOT FAILED_PMRRR)
if(BUILD_TESTS)
  # Build the core tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/core")
  foreach(TEST ${CORE_TESTS})
    add_executable(core-${TEST} ${TEST_DIR}/core/${TEST}.cpp)
    target_link_libraries(core-${TEST} elemental)
    set_target_properties(core-${TEST} 
                          PROPERTIES OUTPUT_NAME ${TEST}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
                          LINK_FLAGS "${BLAS_LAPACK_FLAGS}")
  endforeach(TEST)
  # Build the blas tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/blas")
  foreach(TEST ${BLAS_TESTS})
    add_executable(blas-${TEST} ${TEST_DIR}/blas/${TEST}.cpp)
    target_link_libraries(blas-${TEST} elemental)
    set_target_properties(blas-${TEST} 
                          PROPERTIES OUTPUT_NAME ${TEST}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
                          LINK_FLAGS "${BLAS_LAPACK_FLAGS}")
  endforeach(TEST)
  # Build the lapack tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/lapack")
  foreach(TEST ${LAPACK_TESTS})
    add_executable(lapack-${TEST} ${TEST_DIR}/lapack/${TEST}.cpp)
    target_link_libraries(lapack-${TEST} elemental)
    set_target_properties(lapack-${TEST} 
                          PROPERTIES OUTPUT_NAME ${TEST}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}
                          LINK_FLAGS "${BLAS_LAPACK_FLAGS}")
  endforeach(TEST)
endif(BUILD_TESTS)

