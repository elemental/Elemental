#.rst:
# ElCheckFunctionExists
# ---------------------
#
# Check if a C function can be linked
#
# EL_CHECK_FUNCTION_EXISTS(<function> <variable>)
#
# Check that the <function> is provided by libraries on the system and
# store the result in a <variable>.  This does not verify that any
# system header file declares the function, only that it can be found at
# link time (consider using CheckSymbolExists).
# <variable> will be created as an internal cache variable.
#
# The following variables may be set before calling this macro to modify
# the way the check is run:
#
# ::
#
#   CMAKE_REQUIRED_FLAGS = string of compile command line flags
#   CMAKE_REQUIRED_LINKER_FLAGS = list of linking command line flags
#   CMAKE_REQUIRED_DEFINITIONS = list of macros to define (-DFOO=bar)
#   CMAKE_REQUIRED_INCLUDES = list of include directories
#   CMAKE_REQUIRED_LIBRARIES = list of libraries to link
#   CMAKE_REQUIRED_QUIET = execute quietly without messages

# CMake - Cross Platform Makefile Generator
# Copyright 2000-2015 Kitware, Inc.
# Copyright 2000-2011 Insight Software Consortium
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# 
# * Neither the names of Kitware, Inc., the Insight Software Consortium,
#   nor the names of their contributors may be used to endorse or promote
#   products derived from this software without specific prior written
#   permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
#

# Copyright 2015, Jack Poulson
# All rights reserved.

macro(EL_CHECK_FUNCTION_EXISTS FUNCTION VARIABLE)
  if(NOT DEFINED "${VARIABLE}" OR "x${${VARIABLE}}" STREQUAL "x${VARIABLE}")
    set(MACRO_CHECK_FUNCTION_DEFINITIONS
      "-DCHECK_FUNCTION_EXISTS=${FUNCTION} ${CMAKE_REQUIRED_FLAGS}")
    if(NOT CMAKE_REQUIRED_QUIET)
      message(STATUS "Looking for ${FUNCTION}")
    endif()
    if(CMAKE_REQUIRED_LINKER_FLAGS)
      set(CHECK_FUNCTION_EXISTS_ADD_LINKER_FLAGS -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_REQUIRED_LINKER_FLAGS} -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_REQUIRED_LINKER_FLAGS} -DCMAKE_MODULE_LINKER_FLAGS:STRING=${CMAKE_REQUIRED_LINKER_FLAGS})
    else()
      set(CHECK_FUNCTION_EXISTS_ADD_LINKER_FLAGS)
    endif()
    if(CMAKE_REQUIRED_LIBRARIES)
      set(CHECK_FUNCTION_EXISTS_ADD_LIBRARIES
        LINK_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
    else()
      set(CHECK_FUNCTION_EXISTS_ADD_LIBRARIES)
    endif()
    if(CMAKE_REQUIRED_INCLUDES)
      set(CHECK_FUNCTION_EXISTS_ADD_INCLUDES
        "-DINCLUDE_DIRECTORIES:STRING=${CMAKE_REQUIRED_INCLUDES}")
    else()
      set(CHECK_FUNCTION_EXISTS_ADD_INCLUDES)
    endif()
    try_compile(${VARIABLE}
      ${CMAKE_BINARY_DIR}
      ${CMAKE_ROOT}/Modules/CheckFunctionExists.c
      COMPILE_DEFINITIONS ${CMAKE_REQUIRED_DEFINITIONS}
      ${CHECK_FUNCTION_EXISTS_ADD_LIBRARIES}
      CMAKE_FLAGS -DCOMPILE_DEFINITIONS:STRING=${MACRO_CHECK_FUNCTION_DEFINITIONS} ${CHECK_FUNCTION_EXISTS_ADD_LINKER_FLAGS}
      "${CHECK_FUNCTION_EXISTS_ADD_INCLUDES}"
      OUTPUT_VARIABLE OUTPUT)
    if(${VARIABLE})
      set(${VARIABLE} 1 CACHE INTERNAL "Have function ${FUNCTION}")
      if(NOT CMAKE_REQUIRED_QUIET)
        message(STATUS "Looking for ${FUNCTION} - found")
      endif()
      file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
        "Determining if the function ${FUNCTION} exists passed with the following output:\n"
        "${OUTPUT}\n\n")
    else()
      if(NOT CMAKE_REQUIRED_QUIET)
        message(STATUS "Looking for ${FUNCTION} - not found")
      endif()
      set(${VARIABLE} "" CACHE INTERNAL "Have function ${FUNCTION}")
      file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
        "Determining if the function ${FUNCTION} exists failed with the following output:\n"
        "${OUTPUT}\n\n")
    endif()
  endif()
endmacro()
