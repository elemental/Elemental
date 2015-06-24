#
#  Copyright 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License,
#  which can be found in the LICENSE file in the root directory, or at
#  http://opensource.org/licenses/BSD-2-Clause
#
include(ExternalProject)
include(ElLibraryName)

set(USE_FOUND_PARMETIS FALSE)
if(NOT EL_BUILD_PARMETIS)
  message(STATUS "Searching for previously installed ParMETIS")
  find_package(ParMETIS)
  if(PARMETIS_FOUND)
    include(CheckFunctionExists)
    set(CMAKE_REQUIRED_LIBRARIES ${PARMETIS_LIBRARIES})
    check_function_exists(ParMETIS_ComputeVertexSeparator 
      HAVE_PARMETIS_VERTEXSEP)
    if(HAVE_PARMETIS_VERTEXSEP)
      set(USE_FOUND_PARMETIS TRUE)
    else()
      message(WARNING "ParMETIS was found, but the custom add-on ParMETIS_ComputeVertexSeparator was not, so a custom version of the library must be built")
    endif()
    unset(CMAKE_REQUIRED_LIBRARIES)
  endif()
endif()

if(USE_FOUND_PARMETIS)
  # find_package returns 'PARMETIS_LIBRARIES' but ParMETIS's CMakeLists.txt
  # returns 'PARMETIS_LIBS'
  set(PARMETIS_LIBS ${PARMETIS_LIBRARIES})
else()
  if(NOT DEFINED PARMETIS_URL)
    set(PARMETIS_URL https://github.com/poulson/parmetis.git)
  endif()
  message(STATUS "Will pull ParMETIS from ${PARMETIS_URL}")

  set(PARMETIS_SOURCE_DIR ${PROJECT_BINARY_DIR}/download/parmetis/source)
  set(PARMETIS_BINARY_DIR ${PROJECT_BINARY_DIR}/download/parmetis/build)

  option(METIS_PCRE OFF)
  if(MSVC OR MINGW)
    option(METIS_GKREGEX "Use GKlib's internal regex?" ON)
  else()
    option(METIS_GKREGEX "Use GKlib's internal regex?" OFF)
  endif()

  string(REPLACE ";" "^^" MPI_C_INCSTRING "${MPI_C_INCLUDE_PATH}")
  string(REPLACE ";" "^^" MPI_C_LIBSTRING "${MPI_C_LIBRARIES}")
  ExternalProject_Add(project_parmetis 
    PREFIX ${CMAKE_INSTALL_PREFIX}
    GIT_REPOSITORY ${PARMETIS_URL}
    STAMP_DIR  ${PARMETIS_BINARY_DIR}/stamp
    SOURCE_DIR ${PARMETIS_SOURCE_DIR}
    BINARY_DIR ${PARMETIS_BINARY_DIR}
    TMP_DIR    ${PARMETIS_BINARY_DIR}/tmp
    UPDATE_COMMAND ""
    LIST_SEPARATOR ^^
    CMAKE_ARGS
      -D PCRE=${METIS_PCRE}
      -D GKREGEX=${METIS_GKREGEX}
      -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -D CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
      -D MPI_C_INCLUDE_PATH=${MPI_C_INCLUDE_PATH}
      -D MPI_C_LIBRARIES=${MPI_C_LIBRARIES}
      -D MPI_LINK_FLAGS=${MPI_LINK_FLAGS}
      -D DISABLE_PARMETIS_PROGRAMS=ON
      -D CMAKE_MACOSX_RPATH=${CMAKE_MACOSX_RPATH}
      -D CMAKE_SKIP_BUILD_RPATH=${CMAKE_SKIP_BUILD_RPATH}
      -D CMAKE_BUILD_WITH_INSTALL_RPATH=${CMAKE_BUILD_WITH_INSTALL_RPATH}
      -D CMAKE_INSTALL_RPATH_USE_LINK_PATH=${CMAKE_INSTALL_RPATH_USE_LINK_PATH} 
      -D CMAKE_INSTALL_RPATH=${CMAKE_INSTALL_RPATH}
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
  )

  # Extract the source and install directories
  ExternalProject_Get_Property(project_parmetis source_dir install_dir)
  set(PARMETIS_INCLUDE_DIR "${install_dir}/include")

  # Add targets for libmetis and libparmetis (either shared or static)
  add_library(libmetis ${LIBRARY_TYPE} IMPORTED)
  add_library(libparmetis ${LIBRARY_TYPE} IMPORTED)
  El_library_name(metis_name metis)
  El_library_name(parmetis_name parmetis)
  set(METIS_LIB ${install_dir}/lib/${metis_name})
  set(PARMETIS_LIB ${install_dir}/lib/${parmetis_name})
  set_property(TARGET libmetis PROPERTY IMPORTED_LOCATION ${METIS_LIB})
  set_property(TARGET libparmetis PROPERTY IMPORTED_LOCATION ${PARMETIS_LIB})

  set(PARMETIS_LIBS ${PARMETIS_LIB} ${METIS_LIB})
  set(EL_BUILT_PARMETIS TRUE)
endif()

set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${PARMETIS_LIBS})

message(STATUS "Including ${PARMETIS_INCLUDE_DIR} for external ParMETIS")
include_directories(${PARMETIS_INCLUDE_DIR})
set(EL_HAVE_METIS TRUE)
set(EL_HAVE_PARMETIS TRUE)
