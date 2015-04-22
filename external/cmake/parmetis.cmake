#
#  Copyright 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License,
#  which can be found in the LICENSE file in the root directory, or at
#  http://opensource.org/licenses/BSD-2-Clause
#
include(ExternalProject)

set(USE_FOUND_PARMETIS FALSE)
if(NOT BUILD_PARMETIS)
  find_package(ParMETIS)
  if(PARMETIS_FOUND)
    if(EXISTS "${PARMETIS_DIR}/libparmetis/parmetislib.h")
      set(USE_FOUND_PARMETIS TRUE)
    else()
      message(WARNING "ParMETIS was found, but parmetislib.h was not, and so ParMETIS must be built again to allow for a custom parallel vertex separation routine to be built")
    endif()
  endif()
endif()

if(USE_FOUND_PARMETIS)
  # find_package returns 'PARMETIS_LIBRARIES' but ParMETIS's CMakeLists.txt
  # returns 'PARMETIS_LIBS'
  set(PARMETIS_LIBS ${PARMETIS_LIBRARIES})

  # parmetislib.h is needed for ElParallelBisect, which uses ParMETIS internals
  include_directories(${PARMETIS_DIR}/include ${PARMETIS_DIR}/libparmetis)
else()
  if(NOT DEFINED PARMETIS_URL)
    set(PARMETIS_URL https://github.com/poulson/parmetis.git)
  endif()
  message(STATUS "Will pull ParMETIS from ${PARMETIS_URL}")

  set(PARMETIS_SOURCE_DIR ${PROJECT_BINARY_DIR}/download/parmetis/source)
  set(PARMETIS_BINARY_DIR ${PROJECT_BINARY_DIR}/download/parmetis/build)

  ExternalProject_Add(project_parmetis 
    PREFIX ${CMAKE_INSTALL_PREFIX}
    GIT_REPOSITORY ${PARMETIS_URL}
    STAMP_DIR  ${PARMETIS_BINARY_DIR}/stamp
    SOURCE_DIR ${PARMETIS_SOURCE_DIR}
    BINARY_DIR ${PARMETIS_BINARY_DIR}
    TMP_DIR    ${PARMETIS_BINARY_DIR}/tmp
    CMAKE_ARGS -D CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
               -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
               -D MPI_INCLUDE_PATH=${MPI_C_INCLUDE_PATH}
               -D MPI_LIBRARIES=${MPI_C_LIBRARIES}
               -D MPI_LINK_FLAGS=${MPI_LINK_FLAGS}
               -D DISABLE_PARMETIS_PROGRAMS=ON
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
  )
  add_dependencies(External project_parmetis)

  # Extract the source and install directories
  ExternalProject_Get_Property(project_parmetis source_dir install_dir)

  # Add targets for libmetis and libparmetis (either shared or static)
  if(BUILD_SHARED_LIBS)
    add_library(libmetis SHARED IMPORTED)
    add_library(libparmetis SHARED IMPORTED)
    set_property(TARGET libmetis PROPERTY IMPORTED_LOCATION ${install_dir}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}metis${CMAKE_SHARED_LIBRARY_SUFFIX})
    set_property(TARGET libparmetis PROPERTY IMPORTED_LOCATION ${install_dir}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}parmetis${CMAKE_SHARED_LIBRARY_SUFFIX})
  else()
    add_library(libmetis STATIC IMPORTED)
    add_library(libparmetis STATIC IMPORTED)
    set_property(TARGET libmetis PROPERTY IMPORTED_LOCATION ${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}metis${CMAKE_STATIC_LIBRARY_SUFFIX})
    set_property(TARGET libparmetis PROPERTY IMPORTED_LOCATION ${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}parmetis${CMAKE_STATIC_LIBRARY_SUFFIX})
  endif() 

  # parmetislib.h is needed for ElParallelBisect, which uses ParMETIS internals
  # to construct the vertex separation routine. Furthermore, parmetis includes
  # files from metis/ and metis/GKlib/
  include_directories(${source_dir}/include 
                      ${source_dir}/libparmetis 
                      ${source_dir}/metis/include
                      ${source_dir}/metis/GKlib)

  set(PARMETIS_LIBS libparmetis libmetis)
  set(EL_BUILT_PARMETIS TRUE)
endif()

set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${PARMETIS_LIBS})

set(EL_HAVE_METIS TRUE)
set(EL_HAVE_PARMETIS TRUE)
