#
#  Copyright 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License,
#  which can be found in the LICENSE file in the root directory, or at
#  http://opensource.org/licenses/BSD-2-Clause
#
include(CheckCXXSourceCompiles)

if(OpenMP_C_FLAGS AND OpenMP_CXX_FLAGS)
  set(EL_HAVE_OPENMP TRUE)
  message(STATUS "Using prespecified OpenMP_C_FLAGS=${OpenMP_C_FLAGS}")
  message(STATUS "Using prespecified OpenMP_CXX_FLAGS=${OpenMP_CXX_FLAGS}")
else()
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(EL_HAVE_OPENMP TRUE)
  else()
    set(OpenMP_C_FLAGS "" CACHE STRING "OpenMP C FLAGS")
    set(OpenMP_CXX_FLAGS "" CACHE STRING "OpenMP CXX FLAGS")
    if(EL_HYBRID)
      message(FATAL_ERROR
        "Hybrid build failed because OpenMP support not detected. Please specify OpenMP_C_FLAGS and OpenMP_CXX_FLAGS.")
    endif()
  endif()
endif()

# See if we have 'collapse' support, since Clang is currently setting _OPENMP
# greater than 200805 despite not supporting OpenMP 3.0 (i.e., collapse)
if(EL_HAVE_OPENMP)
  set(CMAKE_REQUIRED_FLAGS ${OpenMP_CXX_FLAGS})
  set(OMP_COLLAPSE_CODE
      "#include <omp.h>
       int main( int argc, char* argv[] ) 
       {
           int k[100];
       #pragma omp collapse(2)
           for( int i=0; i<10; ++i ) 
               for( int j=0; j<10; ++j )
                   k[i+j*10] = i+j; 
           return 0; 
       }")
  check_cxx_source_compiles("${OMP_COLLAPSE_CODE}" EL_HAVE_OMP_COLLAPSE)
  set(CMAKE_REQUIRED_FLAGS)
else()
  set(EL_HAVE_OMP_COLLAPSE FALSE)
endif()
