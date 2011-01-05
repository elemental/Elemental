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
project(PMRRR)
cmake_minimum_required(VERSION 2.8)

# Elemental and PMRRR support four build modes:
#   HybridDebug, HybridRelease, PureDebug, and PureRelease
if(NOT (CMAKE_BUILD_TYPE STREQUAL "HybridDebug") AND
   NOT (CMAKE_BUILD_TYPE STREQUAL "HybridRelease") AND 
   NOT (CMAKE_BUILD_TYPE STREQUAL "PureDebug") AND
   NOT (CMAKE_BUILD_TYPE STREQUAL "PureRelease"))
  message(FATAL_ERROR 
          "CMAKE_BUILD_TYPE must be one of: HybridDebug, HybridRelease, PureDebug, PureRelease")
endif(NOT (CMAKE_BUILD_TYPE STREQUAL "HybridDebug") AND
      NOT (CMAKE_BUILD_TYPE STREQUAL "HybridRelease") AND 
      NOT (CMAKE_BUILD_TYPE STREQUAL "PureDebug") AND
      NOT (CMAKE_BUILD_TYPE STREQUAL "PureRelease"))
set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING 
    "Choose the build type, options are: HybridDebug, HybridRelease, PureDebug, and PureRelease" 
    FORCE)

option(HAVE_SPINLOCKS "Enable if pthread lib supports spinlocks" OFF)
MARK_AS_ADVANCED(HAVE_SPINLOCKS)

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
  set(C_HYBRID_DEBUG_FLAGS "-O2 -g" CACHE STRING
      "C optimization/debug flags for hybrid OpenMP/MPI debug build")
endif(HYBRID_DEBUG)
if(HYBRID_RELEASE)
  set(C_HYBRID_RELEASE_FLAGS "-O3" CACHE STRING
      "C optimization flags for hybrid OpenMP/MPI release build")
endif(HYBRID_RELEASE)
if(PURE_DEBUG)
  set(C_PURE_DEBUG_FLAGS "-O2 -g" CACHE STRING
      "C optimization/debugging flags for pure MPI debug build")
endif(PURE_DEBUG)
if(PURE_RELEASE)
  set(C_PURE_RELEASE_FLAGS "-O3" CACHE STRING
      "C optimization flags for pure MPI release build")
endif(PURE_RELEASE)

# Search the BLAS/LAPACK libs for  to test whether it is a new enough version
# for our purposes. If the LAPACK library is version 3.1 or later, then 
# dstemr will support subset computations. However, dlamch only became 
# threadsafe at version 3.3. If Elemental only requires a pure MPI build, 
# then we don't need to worry about threadsafe routines and can settle for 
# 3.1 or later.
#
# Based upon the netlib changesets, we can determine if the LAPACK version 
# is >= 3.1 by searching for claqr0, and likewise searching for cuncsd to test
# if the version is >= 3.3.
include(CheckFunctionExists)
set(CMAKE_REQUIRED_LIBRARIES ${PURE_BLAS_LAPACK_LIBS})
if(HAVE_DPOTRF)
  check_function_exists(claqr0 LAPACK_3_1)
  check_function_exists(cuncsd LAPACK_3_3)
elseif(HAVE_DPOTRF_POST)
  check_function_exists(claqr0_ LAPACK_3_1)
  check_function_exists(cuncsd_ LAPACK_3_3)
endif(HAVE_DPOTRF)

if(PURE AND NOT LAPACK_3_1)
  set(MISSING_LAPACK_3_1 TRUE)
  message(STATUS "Pure Elemental+PMRRR requires LAPACK >= 3.1.")
else(PURE AND NOT LAPACK_3_1)
  set(MISSING_LAPACK_3_1 FALSE)
endif(PURE AND NOT LAPACK_3_1)

if(HYBRID AND NOT LAPACK_3_3)
  set(MISSING_LAPACK_3_3 TRUE)
  message(STATUS "Hybrid Elemental+PMRRR requires LAPACK >= 3.3.")
elseif(HYBRID AND NOT LAPACK_3_3)
  set(MISSING_LAPACK_3_3 FALSE)
endif(HYBRID AND NOT LAPACK_3_3)

if(MISSING_LAPACK_3_1 OR MISSING_LAPACK_3_3)
  if(NOT CMAKE_Fortran_COMPILER)
    # For some reason enable_language( ... OPTIONAL) causes an error, so we 
    # should instead manually probe for a Fortran compiler first
    #enable_language(Fortran OPTIONAL)
    find_program(FORTRAN_COMPILER NAMES gfortran g77 ifort xlf_r xlf xlf90_r 
                 xlf90 xlf95_r xlf95 xlf2003_r xlf2003)
  endif(NOT CMAKE_Fortran_COMPILER)
  if(CMAKE_Fortran_COMPILER OR FORTRAN_COMPILER)
    enable_language(Fortran)

    if(HYBRID_DEBUG)
      set(Fortran_HYBRID_DEBUG_FLAGS "-O2 -g" CACHE STRING
          "Fortran optimization/debug flags for hybrid threaded/MPI debug build")
    endif(HYBRID_DEBUG)
    if(HYBRID_RELEASE)
      set(Fortran_HYBRID_RELEASE_FLAGS "-O3" CACHE STRING
          "Fortran optimization flags for hybrid threaded/MPI release build")
    endif(HYBRID_RELEASE)
    if(PURE_DEBUG)
      set(Fortran_PURE_DEBUG_FLAGS "-O2 -g" CACHE STRING
          "Fortran optimization/debugging flags for pure MPI debug build")
    endif(PURE_DEBUG)
    if(PURE_RELEASE)
      set(Fortran_PURE_RELEASE_FLAGS "-O3" CACHE STRING
          "Fortran optimization flags for pure MPI release build")
    endif(PURE_RELEASE)

    set(CMAKE_Fortran_FLAGS_HYBRIDDEBUG ${Fortran_HYBRID_DEBUG_FLAGS})
    set(CMAKE_Fortran_FLAGS_HYBRIDRELEASE ${Fortran_HYBRID_RELEASE_FLAGS})
    set(CMAKE_Fortran_FLAGS_PUREDEBUG ${Fortran_PURE_DEBUG_FLAGS})
    set(CMAKE_Fortran_FLAGS_PURERELEASE ${Fortran_PURE_RELEASE_FLAGS})

    # Now determine the way C should interface with our Fortran addons. 
    get_filename_component(Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
    set(LAPACK_ADDONS_LINK_LIBS "")
    if(Fortran_COMPILER_NAME STREQUAL "gfortran")
      if(NOT F_FROM_C_DEFS)
        # gfortran should append an underscore by default
        set(F_FROM_C_DEFS "-DBLAS_POST -DLAPACK_POST")
      endif(NOT F_FROM_C_DEFS)
      if(APPLE)
        string(REGEX REPLACE "bin/gfortran" "lib" GFORTRAN_LIB_DIR
               "${CMAKE_Fortran_COMPILER}")
        find_library(GFORTRAN_LIB gfortran HINTS ${GFORTRAN_LIB_DIR})
        if(GFORTRAN_LIB)
          set(LAPACK_ADDONS_LINK_LIBS ${GFORTRAN_LIB})
        else(GFORTRAN_LIB)
          message(FATAL_ERROR "Could not find gfortran library")
        endif(GFORTRAN_LIB)
      endif(APPLE)
    elseif(Fortran_COMPILER_NAME STREQUAL "bgxlf_r" OR
           Fortran_COMPILER_NAME STREQUAL "bgxlf" OR
           Fortran_COMPILER_NAME STREQUAL "xlf" OR
           Fortran_COMPILER_NAME STREQUAL "xlf_r")
      if(NOT F_FROM_C_DEFS)
        set(F_FROM_C_DEFS "")
      endif(NOT F_FROM_C_DEFS)
    elseif(Fortran_COMPILER_NAME STREQUAL "ifort")
      if(NOT F_FROM_C_DEFS)
        set(F_FROM_C_DEFS "-DBLAS_POST -DLAPACK_POST")
      endif(NOT F_FROM_C_DEFS)
      string(REGEX REPLACE "/bin/([^/]*)/?ifort" "/lib/\\1" IFORT_LIB_DIR
             "${CMAKE_Fortran_COMPILER}")
      find_library(IFCORE_LIB ifcore HINTS ${IFORT_LIB_DIR})
      if(IFCORE_LIB)
        set(LAPACK_ADDONS_LINK_LIBS ${IFCORE_LIB})
      else(IFCORE_LIB)
        message(FATAL_ERROR "Could not find ifort library")
      endif(IFCORE_LIB)
    else(Fortran_COMPILER_NAME STREQUAL "gfortran")
      # The list of other compiler options should grow quickly. It is best      
      # to default to assuming an appended underscore.
      if(NOT F_FROM_C_DEFS) 
        set(F_FROM_C_DEFS "-DBLAS_POST -DLAPACK_POST")
      endif(NOT F_FROM_C_DEFS)
    endif(Fortran_COMPILER_NAME STREQUAL "gfortran")
    set(F_FROM_C_DEFS "${F_FROM_C_DEFS} -DPREPEND_X")
    
    file(GLOB_RECURSE LAPACK_ADDON_SRC RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
         "*.f")
    add_library(lapack-addons STATIC ${LAPACK_ADDON_SRC})
    target_link_libraries(lapack-addons ${LAPACK_ADDONS_LINK_LIBS})
  else(CMAKE_Fortran_COMPILER OR FORTRAN_COMPILER)
    if(MISSING_LAPACK_3_1)
      set(FAILED_PMRRR TRUE)
      set(FAILED_PMRRR TRUE PARENT_SCOPE)
      message(STATUS "Fortran compiler not found; we cannot build PMRRR")
    else(MISSING_LAPACK_3_1)
      set(FAILED_PMRRR TRUE)
      set(FAILED_PMRRR TRUE PARENT_SCOPE)
      message(STATUS "Fortran compiler not found; we cannot build hybrid PMRRR")
    endif(MISSING_LAPACK_3_1)
  endif(CMAKE_Fortran_COMPILER OR FORTRAN_COMPILER)
endif(MISSING_LAPACK_3_1 OR MISSING_LAPACK_3_3)

if(FAILED_PMRRR)
  # Do nothing, we are falling back to building Elemental without PMRRR
else(FAILED_PMRRR)
  # Add in the preprocessor definitions 
  if(MISSING_LAPACK_3_1 OR MISSING_LAPACK_3_3)
    set(C_BASE_FLAGS "${F_FROM_C_DEFS}")
  else(MISSING_LAPACK_3_1 OR MISSING_LAPACK_3_3)
    # Determine the BLAS/LAPACK underscore conventions if necessary
    if(NOT PURE_BLAS_DEFS OR NOT PURE_LAPACK_DEFS)
      set(CMAKE_REQUIRED_LIBRARIES ${PURE_BLAS_LAPACK_LIBS})
      check_function_exists(daxpy HAVE_DAXPY)
      if(HAVE_DAXPY)
        set(PURE_BLAS_DEFS "")
      else(HAVE_DAXPY)
        check_function_exists(daxpy_ HAVE_DAXPY_POST)
        if(HAVE_DAXPY_POST)
          set(PURE_BLAS_DEFS "-DBLAS_POST")
        else(HAVE_DAXPY_POST)
          message(FATAL_ERROR "Could not determine BLAS format.")
        endif(HAVE_DAXPY_POST)
      endif(HAVE_DAXPY)
      check_function_exists(dpotrf HAVE_DPOTRF)
      if(HAVE_DPOTRF)
        set(PURE_LAPACK_DEFS "")
      else(HAVE_DPOTRF)
        check_function_exists(dpotrf_ HAVE_DPOTRF_POST)
        if(HAVE_DPOTRF_POST)
          set(PURE_LAPACK_DEFS "-DLAPACK_POST")
        else(HAVE_DPOTRF_POST)
          message(FATAL_ERROR "Could not determine LAPACK format.")
        endif(HAVE_DPOTRF_POST)
      endif(HAVE_DPOTRF)
    endif(NOT PURE_BLAS_DEFS OR NOT PURE_LAPACK_DEFS)
    set(C_BASE_FLAGS "${PURE_BLAS_DEFS} ${PURE_LAPACK_DEFS}")
  endif(MISSING_LAPACK_3_1 OR MISSING_LAPACK_3_3)
  if(NOT HAVE_SPINLOCKS)
    set(C_BASE_FLAGS "${C_BASE_FLAGS} -DNOSPINLOCKS")
  endif(NOT HAVE_SPINLOCKS)

  # Grab all of the .c and .h files
  file(GLOB_RECURSE PMRRR_C RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.c")
  file(GLOB_RECURSE PMRRR_H RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.h")
  set(PMRRR_SRC "${PMRRR_C};${PMRRR_H}")

  # Try to check if we have an MPI C compiler specified
  if(NOT (CMAKE_C_COMPILER MATCHES "mpi"))
    message(WARNING "C compiler, ${CMAKE_C_COMPILER} does not appear to be an MPI compiler. Please manually set it via the environment variable CC or through the '-DCMAKE_C_COMPILER' cmake flag. You should probably also clean this build folder first.")
  endif(NOT (CMAKE_C_COMPILER MATCHES "mpi"))
    
  # Findpthreads.cmake does not seem to be standard
  #
  #find_package(pthreads REQUIRED)
  #set(C_BASE_FLAGS "${C_BASE_FLAGS} -I${PTHREADS_INCLUDE_DIR}")
  #foreach(DEF ${PTHREADS_DEFINITIONS})
  #  set(C_BASE_FLAGS "${C_BASE_FLAGS} -D${DEF}")
  #endforeach(DEF)
  if(NOT PTHREADS_C_FLAGS)
    # TODO: Make this section more robust by checking a longer list of compile
    #       flags using CMake's builtin routines.
    set(PTHREADS_C_FLAGS "-pthread" CACHE STRING
        "C compilation flags for pthreads support.")
  endif(NOT PTHREADS_C_FLAGS)
  set(C_BASE_FLAGS "${C_BASE_FLAGS} ${PTHREADS_C_FLAGS}")

  # Set the various build flags
  set(CMAKE_C_FLAGS_HYBRIDDEBUG "${C_BASE_FLAGS} ${C_HYBRID_DEBUG_FLAGS}")
  set(CMAKE_C_FLAGS_HYBRIDRELEASE "${C_BASE_FLAGS} ${C_HYBRID_RELEASE_FLAGS}")
  set(CMAKE_C_FLAGS_PUREDEBUG "${C_BASE_FLAGS} ${C_PURE_DEBUG_FLAGS}")
  set(CMAKE_C_FLAGS_PURERELEASE "${C_BASE_FLAGS} ${C_PURE_RELEASE_FLAGS}")

  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
  if(MISSING_LAPACK_3_1 OR MISSING_LAPACK_3_3)
    add_library(pmrrr STATIC ${PMRRR_SRC})
    target_link_libraries(pmrrr lapack-addons)
  else(MISSING_LAPACK_3_1 OR MISSING_LAPACK_3_3)
    add_library(pmrrr STATIC ${PMRRR_SRC})
  endif(MISSING_LAPACK_3_1 OR MISSING_LAPACK_3_3)
endif(FAILED_PMRRR)

