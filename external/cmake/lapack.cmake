#
#  Copyright 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License,
#  which can be found in the LICENSE file in the root directory, or at
#  http://opensource.org/licenses/BSD-2-Clause
#
include(ExternalProject)
include(ElCheckFunctionExists)

if(NOT EL_BUILD_LAPACK)
  find_library(LAPACK NAMES lapack PATHS ${MATH_PATHS})
  if(LAPACK)
    set(CMAKE_REQUIRED_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_REQUIRED_LIBRARIES ${LAPACK} ${MATH_LIBS})
    El_check_function_exists(dstegr  EL_HAVE_DSTEGR)
    El_check_function_exists(dstegr_ EL_HAVE_DSTEGR_POST)
    if(EL_HAVE_DSTEGR)
      set(USE_FOUND_LAPACK TRUE)
      set(EL_LAPACK_SUFFIX)
    elseif(EL_HAVE_DSTEGR_POST)
      set(USE_FOUND_LAPACK TRUE)
      set(EL_LAPACK_SUFFIX _)
    endif() 
    set(CMAKE_REQUIRED_LINKER_FLAGS)
    set(CMAKE_REQUIRED_LIBRARIES)
  endif()
endif()

if(USE_FOUND_LAPACK)
  set(LAPACK_LIBS ${LAPACK})
  set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${LAPACK_LIBS})
  set(EL_HAVE_LAPACK TRUE)
elseif(EL_HAVE_F90_INTERFACE)
  if(NOT DEFINED LAPACK_URL)
    set(LAPACK_URL http://www.netlib.org/lapack/lapack-3.5.0.tgz)
  endif()
  message(STATUS "Will download LAPACK from ${LAPACK_URL}")

  set(LAPACK_SOURCE_DIR ${PROJECT_BINARY_DIR}/download/lapack/source)
  set(LAPACK_BINARY_DIR ${PROJECT_BINARY_DIR}/download/lapack/build)

  ExternalProject_Add(project_lapack
    PREFIX ${CMAKE_INSTALL_PREFIX}
    URL ${LAPACK_URL}
    STAMP_DIR  ${LAPACK_BINARY_DIR}/stamp
    SOURCE_DIR ${LAPACK_SOURCE_DIR}
    BINARY_DIR ${LAPACK_BINARY_DIR}
    TMP_DIR    ${LAPACK_BINARY_DIR}/tmp
    UPDATE_COMMAND ""
    CMAKE_ARGS 
      -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -D CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -D CMAKE_C_FLAGS=${CMAKE_C_FLAGS}
      -D CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
      -D BLAS_LIBRARIES=${MATH_LIBS}
      -D CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
      -D CMAKE_MACOSX_RPATH=${CMAKE_MACOSX_RPATH}
      -D CMAKE_SKIP_BUILD_RPATH=${CMAKE_SKIP_BUILD_RPATH}
      -D CMAKE_BUILD_WITH_INSTALL_RPATH=${CMAKE_BUILD_WITH_INSTALL_RPATH}
      -D CMAKE_INSTALL_RPATH_USE_LINK_PATH=${CMAKE_INSTALL_RPATH_USE_LINK_PATH} 
      -D CMAKE_INSTALL_RPATH=${CMAKE_INSTALL_RPATH}
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
  )
  if(EL_BUILT_BLIS)
    add_dependencies(project_lapack project_blis)
  endif()
  if(EL_BUILT_OPENBLAS)
    add_dependencies(project_lapack project_openblas)
  endif()

  # Extract the source and install directories
  ExternalProject_Get_Property(project_lapack source_dir install_dir)

  # Add a target for liblapack (either shared or static)
  if(BUILD_SHARED_LIBS)
    add_library(liblapack SHARED IMPORTED)
    set(LAPACK_LIB ${install_dir}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}lapack${CMAKE_SHARED_LIBRARY_SUFFIX})
  else()
    add_library(liblapack STATIC IMPORTED)
    set(LAPACK_LIB ${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}lapack${CMAKE_STATIC_LIBRARY_SUFFIX})
  endif() 
  set_property(TARGET liblapack PROPERTY IMPORTED_LOCATION ${LAPACK_LIB})

  set(LAPACK_LIBS ${LAPACK_LIB})
  set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${LAPACK_LIBS})
  set(EL_BUILT_LAPACK TRUE)
  set(EL_HAVE_LAPACK TRUE)
else()
  set(EL_BUILT_LAPACK FALSE)
  set(EL_HAVE_LAPACK FALSE)
endif()
