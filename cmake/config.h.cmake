/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef ELEMENTAL_CONFIG_H
#define ELEMENTAL_CONFIG_H 1

/* Build type and version information */
#define CMAKE_BUILD_TYPE @CMAKE_BUILD_TYPE@
#define Elemental_VERSION_MAJOR @Elemental_VERSION_MAJOR@
#define Elemental_VERSION_MINOR @Elemental_VERSION_MINOR@

/* C compiler info */
#define CMAKE_C_COMPILER    @CMAKE_C_COMPILER@
#define MPI_C_COMPILER      @MPI_C_COMPILER@
#define MPI_C_INCLUDE_PATH  @MPI_C_INCLUDE_PATH@
#define MPI_C_COMPILE_FLAGS @MPI_C_COMPILE_FLAGS@
#define MPI_C_LINK_FLAGS    @MPI_C_LINK_FLAGS@
#define MPI_C_LIBRARIES     @MPI_C_LIBRARIES@

/* C++ compiler info */
#define CMAKE_CXX_COMPILER    @CMAKE_CXX_COMPILER@
#define CXX_FLAGS             @CXX_FLAGS@
#define MPI_CXX_COMPILER      @MPI_CXX_COMPILER@
#define MPI_CXX_INCLUDE_PATH  @MPI_CXX_INCLUDE_PATH@
#define MPI_CXX_COMPILE_FLAGS @MPI_CXX_COMPILE_FLAGS@
#define MPI_CXX_LINK_FLAGS    @MPI_CXX_LINK_FLAGS@
#define MPI_CXX_LIBRARIES     @MPI_CXX_LIBRARIES@

/* Math libraries */
#define MATH_LIBS @MATH_LIBS@

/* Basic configuration options */
#define RESTRICT @RESTRICT@
#cmakedefine RELEASE
#cmakedefine BLAS_POST
#cmakedefine LAPACK_POST
#cmakedefine HAVE_OPENMP
#cmakedefine HAVE_F90_INTERFACE
#cmakedefine WITHOUT_PMRRR
#cmakedefine AVOID_COMPLEX_MPI
#cmakedefine HAVE_FLA_BSVD
#cmakedefine HAVE_REDUCE_SCATTER_BLOCK
#cmakedefine HAVE_MPI_IN_PLACE
#cmakedefine HAVE_MPI3_NONBLOCKING_COLLECTIVES
#cmakedefine HAVE_MPIX_NONBLOCKING_COLLECTIVES
#cmakedefine REDUCE_SCATTER_BLOCK_VIA_ALL_REDUCE
#cmakedefine USE_BYTE_ALLGATHERS

/* Advanced configuration options */
#cmakedefine CACHE_WARNINGS
#cmakedefine UNALIGNED_WARNINGS
#cmakedefine VECTOR_WARNINGS
#cmakedefine POOL_MEMORY
#cmakedefine AVOID_OMP_FMA

#endif // ELEMENTAL_CONFIG_H
