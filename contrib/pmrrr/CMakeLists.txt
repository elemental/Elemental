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

if(CMAKE_BUILD_TYPE)
  message("PMRRR does not use CMAKE_BUILD_TYPE.")
endif(CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE "" CACHE INTERNAL "Removing GUI option" FORCE)

option(HAVE_SPINLOCKS "Enable if pthread lib supports spinlocks" OFF)
MARK_AS_ADVANCED(HAVE_SPINLOCKS)

set(C_OPT_FLAGS "-O3" CACHE STRING "C optimization flags")

set(C_BASE_FLAGS "${BLAS_DEFS} ${LAPACK_DEFS}")
if(NOT HAVE_SPINLOCKS)
  set(C_BASE_FLAGS "${C_BASE_FLAGS} -DNOSPINLOCKS")
endif(NOT HAVE_SPINLOCKS)

# Grab all of the .c and .h files
file(GLOB_RECURSE PMRRR_C RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.c")
file(GLOB_RECURSE PMRRR_H RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.h")
set(PMRRR_SRC "${PMRRR_C};${PMRRR_H}")

# Search for the MPI C compiler
if(MPI_C_COMPILER AND NOT EXISTS ${MPI_C_COMPILER})
  message(FATAL_ERROR 
          "Invalid MPI_C_COMPILER: ${MPI_C_COMPILER} does not exist.")
endif(MPI_C_COMPILER AND NOT EXISTS ${MPI_C_COMPILER})
if(NOT MPI_C_COMPILER)
  find_program(MPI_C_COMPILER mpicc mpixlc mpixlc_r 
    DOC "Full path to MPI C compiler. If blank, CMake will search for it.")
  if(NOT MPI_C_COMPILER)
    message(FATAL_ERROR 
      "Could not find MPI C compiler. Please manually set MPI_C_COMPILER.")
  endif(NOT MPI_C_COMPILER)
endif(NOT MPI_C_COMPILER)
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

# Findpthreads.cmake does not seem to be standard
#
#find_package(pthreads REQUIRED)
#set(C_BASE_FLAGS "${C_BASE_FLAGS} -I${PTHREADS_INCLUDE_DIR}")
#foreach(DEF ${PTHREADS_DEFINITIONS})
#  set(C_BASE_FLAGS "${C_BASE_FLAGS} -D${DEF}")
#endforeach(DEF)
set(C_PTHREADS_FLAGS "-pthread" CACHE STRING 
    "C compile flag for pthreads support.")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
set(CMAKE_C_FLAGS "${C_BASE_FLAGS} ${C_OPT_FLAGS} ${C_PTHREADS_FLAG}")
set(CMAKE_Fortran_FLAGS "${FORTRAN_OPT_FLAGS}")
add_library(pmrrr STATIC ${PMRRR_SRC})

