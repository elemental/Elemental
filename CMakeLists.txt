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
cmake_minimum_required(VERSION 2.8.5)
project(Elemental)

set(Elemental_VERSION_MAJOR 0)
set(Elemental_VERSION_MINOR 71)

# Check if an in-source build was attempted
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "In-source build attempted; please clean the CMake cache and then switch to an out-of-source build, e.g., rm CMakeCache.txt && rm -Rf CMakeFiles/ && mkdir build/ && cd build/ && cmake ..")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Build mode not specified, defaulting to PureRelease build.")
  set(CMAKE_BUILD_TYPE PureRelease)
else(NOT CMAKE_BUILD_TYPE)
  if(NOT (CMAKE_BUILD_TYPE MATCHES HybridDebug) AND
     NOT (CMAKE_BUILD_TYPE MATCHES HybridRelease) AND 
     NOT (CMAKE_BUILD_TYPE MATCHES PureDebug) AND
     NOT (CMAKE_BUILD_TYPE MATCHES PureRelease))
    message(
      FATAL_ERROR 
      "CMAKE_BUILD_TYPE must be in {HybridDebug,HybridRelease,PureDebug,PureRelease}"
    )
  endif(NOT (CMAKE_BUILD_TYPE MATCHES HybridDebug) AND
        NOT (CMAKE_BUILD_TYPE MATCHES HybridRelease) AND 
        NOT (CMAKE_BUILD_TYPE MATCHES PureDebug) AND
        NOT (CMAKE_BUILD_TYPE MATCHES PureRelease))
endif(NOT CMAKE_BUILD_TYPE)

option(BUILD_EXAMPLES "Build the examples" OFF)
option(BUILD_TESTS "Build a collection of test executables" OFF)
if(WIN32)
  option(BUILD_PMRRR "Build and link against eigensolver PMRRR" OFF)
else(WIN32)
  option(BUILD_PMRRR "Build and link against eigensolver PMRRR" ON)
endif(WIN32)

option(WITHOUT_COMPLEX "Avoid all usage of complex numbers" OFF)
if(NOT WITHOUT_COMPLEX)
  option(AVOID_COMPLEX_MPI "Avoid complex MPI routines for robustness" ON)
endif(NOT WITHOUT_COMPLEX)

# Bug fix for Darwin g++ compiler throwing an internal compiler error.
# It seems to have somehow been circumvented after pulling most of the library
# into header files.
option(DISABLE_SCALAR_WRAPPER "Avoid pretty-printing of complex numbers" OFF)

# Bug fix for IBM's compiler miscompiling OpenMP a[i] += alpha b[i] statements
# for complex datatypes
option(AVOID_OMP_FMA "Avoid a bug in the IBM compilers." ON)
mark_as_advanced(AVOID_OMP_FMA)

# Treat the data for MPI_Allgather as unsigned char's in order to avoid known 
# performance bugs
option(USE_BYTE_ALLGATHERS "Avoid BG/P Allgather performance bug." ON)
mark_as_advanced(USE_BYTE_ALLGATHERS)

# This is mainly useful for the library developers
option(CACHE_WARNINGS "Warns when using cache-unfriendly routines" OFF)
mark_as_advanced(CACHE_WARNINGS)

# This is mainly useful for the library developers
option(UNALIGNED_WARNINGS "Warn when performing unaligned redistributions" OFF)
mark_as_advanced(UNALIGNED_WARNINGS)

# This is mainly useful for the library developers
option(VECTOR_WARNINGS "Warn when vector redistribution chances are missed" OFF)
mark_as_advanced(VECTOR_WARNINGS)

# This can easily be performed, but is likely not a good idea
option(POOL_MEMORY "Make Memory class accumulate memory until destruction" OFF)
mark_as_advanced(POOL_MEMORY)

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

if(NOT WIN32)
  if(HYBRID_DEBUG AND NOT CXX_HYBRID_DEBUG_FLAGS)
    set(CXX_HYBRID_DEBUG_FLAGS "-O2 -g" CACHE STRING
        "CXX optimization/debug flags for hybrid OpenMP/MPI debug build")
  endif(HYBRID_DEBUG AND NOT CXX_HYBRID_DEBUG_FLAGS)
  if(HYBRID_RELEASE AND NOT CXX_HYBRID_RELEASE_FLAGS)
    set(CXX_HYBRID_RELEASE_FLAGS "-O3" CACHE STRING
        "CXX optimization flags for hybrid OpenMP/MPI release build")
  endif(HYBRID_RELEASE AND NOT CXX_HYBRID_RELEASE_FLAGS)
  if(PURE_DEBUG AND NOT CXX_PURE_DEBUG_FLAGS)
    set(CXX_PURE_DEBUG_FLAGS "-O2 -g" CACHE STRING
        "CXX optimization/debugging flags for pure MPI debug build")
  endif(PURE_DEBUG AND NOT CXX_PURE_DEBUG_FLAGS)
  if(PURE_RELEASE AND NOT CXX_PURE_RELEASE_FLAGS)
    set(CXX_PURE_RELEASE_FLAGS "-O3" CACHE STRING
        "CXX optimization flags for pure MPI release build")
  endif(PURE_RELEASE AND NOT CXX_PURE_RELEASE_FLAGS)
endif(NOT WIN32)

find_package(MPI)
if(NOT MPI_CXX_FOUND)
  message(FATAL_ERROR "MPI C++ compiler was not found and is required")
endif(NOT MPI_CXX_FOUND)
include_directories(${MPI_CXX_INCLUDE_PATH})

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
        "Hybrid build failed because OpenMP support was not detected. Please try specifying OpenMP_CXX_FLAGS."
      )
    endif(OPENMP_FOUND)
  endif(OpenMP_CXX_FLAGS)
endif(HYBRID)

set(CMAKE_CXX_FLAGS_HYBRIDDEBUG 
    "${CXX_HYBRID_DEBUG_FLAGS} ${OpenMP_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS_HYBRIDRELEASE 
    "${CXX_HYBRID_RELEASE_FLAGS} ${OpenMP_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS_PUREDEBUG 
    "${CXX_PURE_DEBUG_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS_PURERELEASE 
    "${CXX_PURE_RELEASE_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")

# Query the size of a void pointer in order to determine whether or not this is
# a 32-bit system.
# Unfortunately, CMAKE_SIZEOF_VOID_P seems to be incorrect on some Windows 
# platforms, so I am switching to directly checking the size of void*
include(CheckTypeSize)
if(NOT SIZEOF_VOID_P)
  check_type_size("void*" SIZEOF_VOID_P)
  if(${SIZEOF_VOID_P} EQUAL 4)
    set(32_BIT_SYSTEM TRUE)
  else(${SIZEOF_VOID_P} EQUAL 4)
    set(32_BIT_SYSTEM FALSE)    
  endif(${SIZEOF_VOID_P} EQUAL 4)
endif(NOT SIZEOF_VOID_P)

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
      set(INTEL_ARCH_SUBDIRS 32/lib lib/32)
      set(MKL_ARCH_SUBDIRS 32/lib lib/32)
    else(32_BIT_SYSTEM)
      set(INTEL_ARCH_SUBDIRS em64t/lib intel64/lib lib/em64t lib/intel64)
      set(MKL_ARCH_SUBDIRS em64t/lib intel64/lib lib/em64t lib/intel64 
          mkl/em64t/lib mkl/intel64/lib mkl/lib/em64t mkl/lib/intel64)
    endif(32_BIT_SYSTEM)
    if(INTEL_ROOT)
      message(STATUS "Will also search INTEL_ROOT=${INTEL_ROOT} for MKL")
      message(STATUS "INTEL_ARCH_SUBDIRS=${INTEL_ARCH_SUBDIRS}")
    endif(INTEL_ROOT)
    if(MKL_ROOT)
      message(STATUS "Will also search MKL_ROOT=${MKL_ROOT} for MKL")
      message(STATUS "MKL_ARCH_SUBDIRS=${INTEL_ARCH_SUBDIRS}")
    endif(MKL_ROOT)
    # There is a cyclic dependency between MKL_CORE and 
    # MKL_SEQUENTIAL/MKL_INTEL_THREAD with some versions of MKL
    if(NOT WIN32)
      find_library(PTHREAD_LIB pthread DOC "Pthreads library")
    endif(NOT WIN32)
    find_library(
      CMATH_LIB
      NAMES m libmmd.lib
      PATHS ${INTEL_ROOT} 
      PATH_SUFFIXES ${INTEL_ARCH_SUBDIRS}
      DOC "C math library"
    )
    find_library(
      IOMP5_LIB
      NAMES libiomp5.a libiomp5md.lib
      PATHS ${INTEL_ROOT} 
      PATH_SUFFIXES ${INTEL_ARCH_SUBDIRS}
    )
    find_library(
      GUIDE_LIB
      NAMES libguide.a libguide.lib
      PATHS ${INTEL_ROOT}
      PATH_SUFFIXES ${INTEL_ARCH_SUBDIRS}
    )
    find_library(
      MKL_CORE_LIB
      NAMES libmkl_core.a mkl_core.lib 
      PATHS ${MKL_ROOT}
      PATH_SUFFIXES ${MKL_ARCH_SUBDIRS}
    )
    find_library(
      MKL_INTEL_LIB
      NAMES libmkl_intel_lp64.a mkl_intel_lp64.lib 
      PATHS ${MKL_ROOT}
      PATH_SUFFIXES ${MKL_ARCH_SUBDIRS}
    )
    if(PURE)
      find_library(
        MKL_SEQUENTIAL_LIB
        NAMES libmkl_sequential.a mkl_sequential.lib
        PATHS ${MKL_ROOT}
        PATH_SUFFIXES ${MKL_ARCH_SUBDIRS}
      )
      set(MKL_REQUIRED MKL_INTEL MKL_SEQUENTIAL MKL_CORE MKL_SEQUENTIAL CMATH
                       GUIDE IOMP5)
    else(PURE)
      find_library(
        MKL_INTEL_THREAD_LIB
        NAMES libmkl_intel_thread.a mkl_intel_thread.lib
        PATHS ${MKL_ROOT}
        PATH_SUFFIXES ${MKL_ARCH_SUBDIRS}
      )
      set(MKL_REQUIRED MKL_INTEL MKL_INTEL_THREAD MKL_CORE MKL_INTEL_THREAD 
                       CMATH GUIDE IOMP5)
    endif(PURE)
    if(NOT WIN32)
      list(APPEND MKL_REQUIRED PTHREAD)
    endif(NOT WIN32)
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
      message(STATUS "MATH_LIBS=${MATH_LIBS}")
      set(NEED_MATH FALSE)
    endif(MKL_FOUND)
    if(NEED_MATH)
      # Look for default BLAS and LAPACK
      if(REFERENCE_ROOT)
        message(STATUS 
          "Will search REFERENCE_ROOT=${REFERENCE_ROOT} for math libs")
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
        message(WARNING "Using reference BLAS/LAPACK.")
      else(REFERENCE_FOUND) 
        message(FATAL_ERROR 
          "Could not find BLAS/LAPACK libraries. Please specify MATH_LIBS")
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

# Look for restrict support
include(CheckCXXSourceCompiles)
set(RESTRICT_CODE
    "int main(void)
     {
         int* RESTRICT a;
         return 0;
     }
    ")
set(CMAKE_REQUIRED_DEFINITIONS "-DRESTRICT=__restrict__")
check_cxx_source_compiles("${RESTRICT_CODE}" HAVE___restrict__)
if(HAVE___restrict__)
  set(RESTRICT "__restrict__")
  message(STATUS "Using __restrict__ keyword.")
else(HAVE___restrict__)
  set(CMAKE_REQUIRED_DEFINITIONS "-DRESTRICT=__restrict")
  check_cxx_source_compiles("${RESTRICT_CODE}" HAVE___restrict)
  if(HAVE___restrict)
    set(RESTRICT "__restrict")
    message(STATUS "Using __restrict keyword.")
  else(HAVE___restrict)
    set(CMAKE_REQUIRED_DEFINITIONS "-DRESTRICT=restrict")
    check_cxx_source_compiles("${RESTRICT_CODE}" HAVE_restrict)
    if(HAVE_restrict)
      set(RESTRICT "restrict")
      message(STATUS "Using restrict keyword.")
    else(HAVE_restrict)
      set(RESTRICT "")
      message(STATUS "Could not find a restrict keyword.")
    endif(HAVE_restrict)
  endif(HAVE___restrict)
endif(HAVE___restrict__)

add_subdirectory(external/plcg)
if(BUILD_PMRRR)
  add_subdirectory(external/pmrrr)
endif(BUILD_PMRRR)

if(BUILD_PMRRR AND NOT FAILED_PMRRR)
  set(WITHOUT_PMRRR FALSE)
else(BUILD_PMRRR AND NOT FAILED_PMRRR)
  set(WITHOUT_PMRRR TRUE)
endif(BUILD_PMRRR AND NOT FAILED_PMRRR)

# Create the Elemental configuration header
configure_file( 
  ${PROJECT_SOURCE_DIR}/include/elemental/config.h.cmake
  ${PROJECT_BINARY_DIR}/include/elemental/config.h 
) 
install(FILES ${PROJECT_BINARY_DIR}/include/elemental/config.h 
        DESTINATION include/elemental)

# Grab all of the .c, .cpp, .h, and .hpp Elemental files
file(GLOB_RECURSE ELEMENTAL_CPP RELATIVE ${PROJECT_SOURCE_DIR} 
     "src/*.c" "src/*.cpp")
file(GLOB_RECURSE ELEMENTAL_HEADERS RELATIVE ${PROJECT_SOURCE_DIR} 
     "include/*.h" "include/*.hpp")
set(ELEMENTAL_SRC "${ELEMENTAL_CPP};${ELEMENTAL_HEADERS}")

# Create a dummy library in order to be able to force the math (and potentially
# MPI libraries) to be linked last
add_library(cmake-dummy-lib STATIC cmake/CMakeDummyFunction.cpp)
target_link_libraries(cmake-dummy-lib ${MATH_LIBS} ${MPI_CXX_LIBRARIES})

# The main library
if(BUILD_PMRRR AND NOT FAILED_PMRRR)
  add_library(elemental STATIC ${ELEMENTAL_SRC})
  target_link_libraries(elemental pmrrr plcg)
  if(MPI_CXX_LINK_FLAGS AND MPI_C_LINK_FLAGS)
    set(MPI_LINK_FLAGS "${MPI_CXX_LINK_FLAGS} ${MPI_C_LINK_FLAGS}")
  elseif(MPI_CXX_LINK_FLAGS)
    set(MPI_LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  elseif(MPI_C_LINK_FLAGS)
    set(MPI_LINK_FLAGS "${MPI_C_LINK_FLAGS}")
  endif(MPI_CXX_LINK_FLAGS AND MPI_C_LINK_FLAGS)
else(BUILD_PMRRR AND NOT FAILED_PMRRR)
  add_library(elemental STATIC ${ELEMENTAL_SRC})
  target_link_libraries(elemental plcg)
  set(MPI_LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
endif(BUILD_PMRRR AND NOT FAILED_PMRRR)
install(TARGETS elemental DESTINATION lib)

# Define the header-file preparation rules
set(PREPARED_HEADERS)
foreach(HEADER ${ELEMENTAL_HEADERS})
  add_custom_command(
    OUTPUT ${PROJECT_BINARY_DIR}/${HEADER}
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/${HEADER}
            ${PROJECT_BINARY_DIR}/${HEADER}
    DEPENDS "${PROJECT_SOURCE_DIR}/${HEADER}"
  )
  list(APPEND PREPARED_HEADERS ${PROJECT_BINARY_DIR}/${HEADER})

  get_filename_component(HEADER_PATH ${HEADER} PATH)
  install(FILES ${PROJECT_BINARY_DIR}/${HEADER} DESTINATION ${HEADER_PATH})
endforeach(HEADER)
add_custom_target(prepare_elemental_headers DEPENDS ${PREPARED_HEADERS})
add_dependencies(elemental prepare_elemental_headers)

# Make sure the Elemental headers can be found
include_directories("${PROJECT_BINARY_DIR}/include")

# Build the test drivers if necessary
if(BUILD_TESTS)
  set(TEST_DIR ${PROJECT_SOURCE_DIR}/tests)

  set(CORE_TESTS AxpyInterface DifferentGrids DistMatrix Matrix)
  set(BASIC_TESTS Gemm Hemm Her2k Herk Symm Symv Syr2k Syrk Trmm Trsm)
  set(ADVANCED_TESTS ApplyPackedReflectors Cholesky Hegst HermitianTridiag LDL 
                     LU LQ QR TriangularInverse)

  if(BUILD_PMRRR AND NOT FAILED_PMRRR)
    list(APPEND ADVANCED_TESTS HermitianGenDefiniteEig)
    list(APPEND ADVANCED_TESTS HermitianEig)
  endif(BUILD_PMRRR AND NOT FAILED_PMRRR)

  # Build the core tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/tests/core")
  foreach(TEST ${CORE_TESTS})
    add_executable(tests-core-${TEST} ${TEST_DIR}/core/${TEST}.cpp)
    target_link_libraries(tests-core-${TEST} elemental cmake-dummy-lib)
    set_target_properties(tests-core-${TEST} PROPERTIES
                          OUTPUT_NAME ${TEST}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(tests-core-${TEST} PROPERTIES
                            LINK_FLAGS ${MPI_LINK_FLAGS})
    endif(MPI_LINK_FLAGS)
    install(TARGETS tests-core-${TEST} DESTINATION bin/tests/core)
  endforeach(TEST)
  # Build the basic tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/tests/basic")
  foreach(TEST ${BASIC_TESTS})
    add_executable(tests-basic-${TEST} ${TEST_DIR}/basic/${TEST}.cpp)
    target_link_libraries(tests-basic-${TEST} elemental cmake-dummy-lib)
    set_target_properties(tests-basic-${TEST} PROPERTIES
                          OUTPUT_NAME ${TEST}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(tests-basic-${TEST} PROPERTIES
                            LINK_FLAGS ${MPI_LINK_FLAGS})
    endif(MPI_LINK_FLAGS)
    install(TARGETS tests-basic-${TEST} DESTINATION bin/tests/basic)
  endforeach(TEST)
  # Build the advanced tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/tests/advanced")
  foreach(TEST ${ADVANCED_TESTS})
    add_executable(tests-advanced-${TEST} ${TEST_DIR}/advanced/${TEST}.cpp)
    target_link_libraries(tests-advanced-${TEST} elemental cmake-dummy-lib)
    set_target_properties(tests-advanced-${TEST} PROPERTIES
                          OUTPUT_NAME ${TEST}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(tests-advanced-${TEST} PROPERTIES
                            LINK_FLAGS ${MPI_LINK_FLAGS})
    endif(MPI_LINK_FLAGS)
    install(TARGETS tests-advanced-${TEST} DESTINATION bin/tests/advanced)
  endforeach(TEST)
  # Build the export tests
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/tests/exports/c")
  foreach(TEST ${EXPORT_TESTS})
    add_executable(tests-exports-${TEST} ${TEST_DIR}/exports/${TEST}.c)
    target_link_libraries(tests-exports-${TEST} elemental cmake-dummy-lib)
    set_target_properties(tests-exports-${TEST} PROPERTIES
                          OUTPUT_NAME ${TEST}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(tests-exports-${TEST} PROPERTIES
                            LINK_FLAGS ${MPI_LINK_FLAGS})
    endif(MPI_LINK_FLAGS)
    install(TARGETS tests-exports-c-${TEST} DESTINATION bin/tests/exports/c)
  endforeach(TEST)
endif(BUILD_TESTS)

# Build the example drivers if necessary
if(BUILD_EXAMPLES)
  set(EXAMPLE_DIR ${PROJECT_SOURCE_DIR}/examples)

  set(CORE_EXAMPLES Identity Random)
  set(ADVANCED_EXAMPLES)

  if(BUILD_PMRRR AND NOT FAILED_PMRRR)
    if(NOT WITHOUT_COMPLEX)
      list(APPEND ADVANCED_EXAMPLES HermitianEig 
                                    ComplexHermitianFunction 
                                    RealHermitianFunction
                                    HermitianPseudoinverse )
    endif(NOT WITHOUT_COMPLEX)
    list(APPEND ADVANCED_EXAMPLES SymmetricEig RealSymmetricFunction)
  endif(BUILD_PMRRR AND NOT FAILED_PMRRR)

  # Build the core test(s)
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/examples/core")
  foreach(EXAMPLE ${CORE_EXAMPLES})
    add_executable(examples-core-${EXAMPLE} ${EXAMPLE_DIR}/core/${EXAMPLE}.cpp)
    target_link_libraries(examples-core-${EXAMPLE} elemental cmake-dummy-lib)
    set_target_properties(examples-core-${EXAMPLE} PROPERTIES 
                          OUTPUT_NAME ${EXAMPLE} 
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(examples-core-${EXAMPLE} PROPERTIES
                            LINK_FLAGS ${MPI_LINK_FLAGS})
    endif(MPI_LINK_FLAGS)
    install(TARGETS examples-core-${EXAMPLE} DESTINATION bin/examples/core)
  endforeach(EXAMPLE)

  # Build the advanced test(s)
  set(OUTPUT_DIR "${PROJECT_BINARY_DIR}/bin/examples/advanced")
  foreach(EXAMPLE ${ADVANCED_EXAMPLES})
    add_executable(examples-advanced-${EXAMPLE} 
                   ${EXAMPLE_DIR}/advanced/${EXAMPLE}.cpp)
    target_link_libraries(examples-advanced-${EXAMPLE} 
                          elemental cmake-dummy-lib)
    set_target_properties(examples-advanced-${EXAMPLE} PROPERTIES
                          OUTPUT_NAME ${EXAMPLE}
                          RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR})
    if(MPI_LINK_FLAGS)
      set_target_properties(examples-advanced-${EXAMPLE} PROPERTIES
                            LINK_FLAGS ${MPI_LINK_FLAGS})
    endif(MPI_LINK_FLAGS)
    install(TARGETS examples-advanced-${EXAMPLE} 
            DESTINATION bin/examples/advanced)
  endforeach(EXAMPLE)
endif(BUILD_EXAMPLES)

# Simplify the inclusion of Elemental as a subproject in another build system
#
#set(MPI_C_COMPILER ${MPI_C_COMPILER} PARENT_SCOPE)
#set(MPI_C_INCLUDE_PATH ${MPI_C_INCLUDE_PATH} PARENT_SCOPE)
#set(MPI_C_COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS} PARENT_SCOPE)
#set(MPI_C_LINK_FLAGS ${MPI_C_LINK_FLAGS} PARENT_SCOPE)
#set(MPI_C_LIBRARIES ${MPI_C_LIBRARIES} PARENT_SCOPE)
#set(MPI_CXX_COMPILER ${MPI_CXX_COMPILER} PARENT_SCOPE)
#set(MPI_CXX_INCLUDE_PATH ${MPI_CXX_INCLUDE_PATH} PARENT_SCOPE)
#set(MPI_CXX_COMPILE_FLAGS ${MPI_CXX_COMPILE_FLAGS} PARENT_SCOPE)
#set(MPI_CXX_LINK_FLAGS ${MPI_CXX_LINK_FLAGS} PARENT_SCOPE)
#set(MPI_CXX_LIBRARIES ${MPI_CXX_LIBRARIES} PARENT_SCOPE)
#set(MATH_LIBS ${MATH_LIBS} PARENT_SCOPE)
#set(RESTRICT ${RESTRICT} PARENT_SCOPE)
#set(RELEASE ${RELEASE} PARENT_SCOPE)
#set(BLAS_POST ${BLAS_POST} PARENT_SCOPE)
#set(LAPACK_POST ${LAPACK_POST} PARENT_SCOPE)
#set(WITHOUT_PMRRR ${WITHOUT_PMRRR} PARENT_SCOPE)
#set(DISABLE_SCALAR_WRAPPER ${DISABLE_SCALAR_WRAPPER} PARENT_SCOPE)
#set(WITHOUT_COMPLEX ${WITHOUT_COMPLEX} PARENT_SCOPE)
#set(AVOID_COMPLEX_MPI ${AVOID_COMPLEX_MPI} PARENT_SCOPE)
#set(USE_BYTE_ALLGATHERS ${USE_BYTE_ALLGATHERS} PARENT_SCOPE)
