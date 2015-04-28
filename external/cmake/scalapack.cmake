#
#  Copyright 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License,
#  which can be found in the LICENSE file in the root directory, or at
#  http://opensource.org/licenses/BSD-2-Clause
#
include(ExternalProject)

set(USE_FOUND_SCALAPACK FALSE)
if(NOT EL_BUILD_SCALAPACK)
  find_library(ScaLAPACK NAMES scalapack PATHS ${MATH_PATHS})
  if(ScaLAPACK)
    get_filename_component(ScaLAPACK_DIR ScaLAPACK DIRECTORY)
    if(EXISTS "${ScaLAPACK_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}blacs${CMAKE_SHARED_LIBRARY_SUFFIX}" OR 
       EXISTS "${ScaLAPACK_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}blacs${CMAKE_STATIC_LIBRARY_SUFFIX}")
      message(WARNING "Found blacs library in ${ScaLAPACK_DIR}; ScaLAPACK appears to be too old to link to")
    else()
      set(CMAKE_REQUIRED_FLAGS "${MPI_C_COMPILE_FLAGS} ${MPI_LINK_FLAGS}")
      set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
      set(CMAKE_REQUIRED_LIBRARIES ${ScaLAPACK} ${MATH_LIBS} ${MPI_C_LIBRARIES})
      # NOTE: pdsyngst was chosen because MKL ScaLAPACK only defines pdsyngst_,
      #       but not pdsyngst, despite defining both pdpotrf and pdpotrf_. 
      check_function_exists(pdsyngst  EL_HAVE_PDSYNGST)
      check_function_exists(pdsyngst_ EL_HAVE_PDSYNGST_POST)
      check_function_exists(Csys2blacs_handle EL_HAVE_CSYS2BLACS)
      if(EL_HAVE_PDSYNGST)
        check_function_exists(pdlaqr0 EL_HAVE_PDLAQR0)
        check_function_exists(pdlaqr1 EL_HAVE_PDLAQR1)
        if(NOT EL_HAVE_PDLAQR0 OR NOT EL_HAVE_PDLAQR1 OR 
           NOT EL_HAVE_CSYS2BLACS)
          message(STATUS "Found ${ScaLAPACK}, but PDLAQR{0,1} and Csys2blacs_handle were not supported.")
        else()
          set(USE_FOUND_SCALAPACK TRUE)
          set(EL_HAVE_SCALAPACK_SUFFIX FALSE)
        endif()
      elseif(EL_HAVE_PDSYNGST_POST)
        check_function_exists(pdlaqr0_ EL_HAVE_PDLAQR0_POST)
        check_function_exists(pdlaqr1_ EL_HAVE_PDLAQR1_POST)
        if(NOT EL_HAVE_PDLAQR0_POST OR NOT EL_HAVE_PDLAQR1_POST OR
           NOT EL_HAVE_CSYS2BLACS)
          message(STATUS "Found ${ScaLAPACK} but PDLAQR{0,1} and Csys2blacs_handle were not supported.")
        else()
          set(USE_FOUND_SCALAPACK TRUE)
          set(EL_SCALAPACK_SUFFIX _)
        endif()
      endif() 
    endif()
  endif()
endif()

if(USE_FOUND_SCALAPACK)
  set(SCALAPACK_LIBS ${ScaLAPACK})
  set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${SCALAPACK_LIBS})
  set(EL_HAVE_SCALAPACK TRUE)
elseif(EL_HAVE_F90_INTERFACE)
  if(NOT DEFINED SCALAPACK_URL)
    set(SCALAPACK_URL https://github.com/poulson/scalapack.git)
  endif()
  message(STATUS "Will pull ScaLAPACK from ${SCALAPACK_URL}")

  set(SCALAPACK_SOURCE_DIR ${PROJECT_BINARY_DIR}/download/scalapack/source)
  set(SCALAPACK_BINARY_DIR ${PROJECT_BINARY_DIR}/download/scalapack/build)

  ExternalProject_Add(project_scalapack
    PREFIX ${CMAKE_INSTALL_PREFIX}
    GIT_REPOSITORY ${SCALAPACK_URL}
    STAMP_DIR  ${SCALAPACK_BINARY_DIR}/stamp
    SOURCE_DIR ${SCALAPACK_SOURCE_DIR}
    BINARY_DIR ${SCALAPACK_BINARY_DIR}
    TMP_DIR    ${SCALAPACK_BINARY_DIR}/tmp
    CMAKE_ARGS 
      -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -D CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -D CMAKE_C_FLAGS=${CMAKE_C_FLAGS}
      -D CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
      -D MPI_C_INCLUDE_PATH=${MPI_C_INCLUDE_PATH}
      -D MPI_Fortran_INCLUDE_PATH=${MPI_Fortran_INCLUDE_PATH}
      -D MPI_C_COMPILE_FLAGS=${MPI_C_COMPILE_FLAGS}
      -D MPI_Fortran_COMPILE_FLAGS=${MPI_Fortran_COMPILE_FLAGS}
      -D MPI_C_LIBRARIES=${MPI_C_LIBRARIES}
      -D MPI_Fortran_LIBRARIES=${MPI_Fortran_LIBRARIES}
      -D MPI_LINK_FLAGS=${MPI_LINK_FLAGS}
      -D LAPACK_LIBRARIES=${MATH_LIBS}
      -D CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
      -D CMAKE_MACOSX_RPATH=${CMAKE_MACOSX_RPATH}
      -D CMAKE_SKIP_BUILD_RPATH=${CMAKE_SKIP_BUILD_RPATH}
      -D CMAKE_BUILD_WITH_INSTALL_RPATH=${CMAKE_BUILD_WITH_INSTALL_RPATH}
      -D CMAKE_INSTALL_RPATH_USE_LINK_PATH=${CMAKE_INSTALL_RPATH_USE_LINK_PATH} 
      -D CMAKE_INSTALL_RPATH=${CMAKE_INSTALL_RPATH}
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
  )
  add_dependencies(External project_scalapack)

  # Extract the source and install directories
  ExternalProject_Get_Property(project_scalapack source_dir install_dir)

  # Add targets for libmetis and libparmetis (either shared or static)
  if(BUILD_SHARED_LIBS)
    add_library(libscalapack SHARED IMPORTED)
    set_property(TARGET libscalapack PROPERTY IMPORTED_LOCATION ${install_dir}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}scalapack${CMAKE_SHARED_LIBRARY_SUFFIX})
  else()
    add_library(libscalapack STATIC IMPORTED)
    set_property(TARGET libscalapack PROPERTY IMPORTED_LOCATION ${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}scalapack${CMAKE_STATIC_LIBRARY_SUFFIX})
  endif() 

  set(SCALAPACK_LIBS libscalapack)
  set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${SCALAPACK_LIBS})
  set(EL_BUILT_SCALAPACK TRUE)
  set(EL_HAVE_SCALAPACK TRUE)
else()

endif()


