/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef ELEM_CONFIG_H
#define ELEM_CONFIG_H

/* Build type and version information */
#define ELEM_GIT_SHA1 "@GIT_SHA1@"
#define Elemental_VERSION_MAJOR "@Elemental_VERSION_MAJOR@"
#define Elemental_VERSION_MINOR "@Elemental_VERSION_MINOR@"
#define ELEM_CMAKE_BUILD_TYPE "@CMAKE_BUILD_TYPE@"
#cmakedefine ELEM_RELEASE

/* C compiler info */
#define ELEM_CMAKE_C_COMPILER    "@CMAKE_C_COMPILER@"
#define ELEM_MPI_C_COMPILER      "@MPI_C_COMPILER@"
#define ELEM_MPI_C_INCLUDE_PATH  "@MPI_C_INCLUDE_PATH@"
#define ELEM_MPI_C_COMPILE_FLAGS "@MPI_C_COMPILE_FLAGS@"
#define ELEM_MPI_C_LINK_FLAGS    "@MPI_C_LINK_FLAGS@"
#define ELEM_MPI_C_LIBRARIES     "@MPI_C_LIBRARIES@"

/* C++ compiler info */
#define ELEM_CMAKE_CXX_COMPILER    "@CMAKE_CXX_COMPILER@"
#define ELEM_CXX_FLAGS             "@CXX_FLAGS@"
#define ELEM_MPI_CXX_COMPILER      "@MPI_CXX_COMPILER@"
#define ELEM_MPI_CXX_INCLUDE_PATH  "@MPI_CXX_INCLUDE_PATH@"
#define ELEM_MPI_CXX_COMPILE_FLAGS "@MPI_CXX_COMPILE_FLAGS@"
#define ELEM_MPI_CXX_LINK_FLAGS    "@MPI_CXX_LINK_FLAGS@"
#define ELEM_MPI_CXX_LIBRARIES     "@MPI_CXX_LIBRARIES@"

/* Math libraries */
#define ELEM_MATH_LIBS "@MATH_LIBS@"
#cmakedefine ELEM_BLAS_POST
#cmakedefine ELEM_LAPACK_POST
#cmakedefine ELEM_HAVE_SCALAPACK
#cmakedefine ELEM_SCALAPACK_POST
#cmakedefine ELEM_HAVE_FLA_BSVD
#define ELEM_FORT_LOGICAL @ELEM_FORT_LOGICAL@
#define ELEM_FORT_TRUE    @ELEM_FORT_TRUE@
#define ELEM_FORT_FALSE   @ELEM_FORT_FALSE@

/* Basic configuration options */
#define ELEM_RESTRICT @RESTRICT@
#cmakedefine ELEM_HAVE_OPENMP
#cmakedefine ELEM_HAVE_OMP_COLLAPSE
#cmakedefine ELEM_HAVE_QT5
#cmakedefine ELEM_HAVE_F90_INTERFACE
#cmakedefine ELEM_AVOID_COMPLEX_MPI
#cmakedefine ELEM_HAVE_CXX11RANDOM
#cmakedefine ELEM_HAVE_STEADYCLOCK
#cmakedefine ELEM_HAVE_NOEXCEPT
#cmakedefine ELEM_HAVE_MPI_REDUCE_SCATTER_BLOCK
#cmakedefine ELEM_HAVE_MPI_IN_PLACE
#cmakedefine ELEM_HAVE_MPI_LONG_LONG
#cmakedefine ELEM_HAVE_MPI_COMM_SET_ERRHANDLER
#cmakedefine ELEM_HAVE_MPI_INIT_THREAD
#cmakedefine ELEM_HAVE_MPI_QUERY_THREAD
#cmakedefine ELEM_HAVE_MPI3_NONBLOCKING_COLLECTIVES
#cmakedefine ELEM_HAVE_MPIX_NONBLOCKING_COLLECTIVES
#cmakedefine ELEM_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
#cmakedefine ELEM_USE_BYTE_ALLGATHERS
#cmakedefine ELEM_USE_64BIT_INTS

/* Advanced configuration options */
#cmakedefine ELEM_ZERO_INIT
#cmakedefine ELEM_CACHE_WARNINGS
#cmakedefine ELEM_UNALIGNED_WARNINGS
#cmakedefine ELEM_VECTOR_WARNINGS
#cmakedefine ELEM_POOL_MEMORY
#cmakedefine ELEM_AVOID_OMP_FMA

#cmakedefine ELEM_HAVE_VALGRIND

#endif /* ELEM_CONFIG_H */
