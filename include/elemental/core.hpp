/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_HPP
#define ELEM_CORE_HPP

#include "mpi.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <random>
#include <vector>

// The DEBUG_ONLY macro is, to the best of my knowledge, the only preprocessor
// name defined by Elemental that is not namespaced with "ELEM". Given how
// frequently it is used, I will leave it as-is unless/until a user/developer 
// complains.
#ifdef ELEM_RELEASE
# define DEBUG_ONLY(cmd) 
#else
# define DEBUG_ONLY(cmd) cmd;
#endif

// If defined, the _OPENMP macro contains the date of the specification
#ifdef ELEM_HAVE_OPENMP
# include <omp.h>
# define ELEM_PARALLEL_FOR _Pragma("omp parallel for")
# ifdef ELEM_HAVE_OMP_COLLAPSE
#  define ELEM_PARALLEL_FOR_COLLAPSE2 _Pragma("omp parallel for collapse(2)")
# else
#  define ELEM_PARALLEL_FOR_COLLAPSE2 ELEM_PARALLEL_FOR
# endif
#else
# define ELEM_PARALLEL_FOR 
# define ELEM_PARALLEL_FOR_COLLAPSE2
#endif

#ifdef ELEM_AVOID_OMP_FMA
# define ELEM_FMA_PARALLEL_FOR 
#else
# define ELEM_FMA_PARALLEL_FOR ELEM_PARALLEL_FOR
#endif
#ifdef ELEM_PARALLELIZE_INNER_LOOPS
# define ELEM_INNER_PARALLEL_FOR           ELEM_PARALLEL_FOR
# define ELEM_INNER_PARALLEL_FOR_COLLAPSE2 ELEM_PARALLEL_FOR_COLLAPSE2
# define ELEM_OUTER_PARALLEL_FOR 
# define ELEM_OUTER_PARALLEL_FOR_COLLAPSE2
#else
# define ELEM_INNER_PARALLEL_FOR
# define ELEM_INNER_PARALLEL_FOR_COLLAPSE2
# define ELEM_OUTER_PARALLEL_FOR           ELEM_PARALLEL_FOR
# define ELEM_OUTER_PARALLEL_FOR_COLLAPSE2 ELEM_PARALLEL_FOR_COLLAPSE2
#endif

#ifdef ELEM_HAVE_NOEXCEPT
# define ELEM_NOEXCEPT noexcept
#else
# define ELEM_NOEXCEPT
#endif

#if defined(ELEM_BLAS_POST)
# define ELEM_BLAS(name) name ## _
#else
# define ELEM_BLAS(name) name
#endif

#if defined(ELEM_LAPACK_POST)
# define ELEM_LAPACK(name) name ## _
#else
# define ELEM_LAPACK(name) name
#endif

#if defined(ELEM_HAVE_SCALAPACK)
# if defined(ELEM_SCALAPACK_POST)
#  define ELEM_SCALAPACK(name) name ## _
# else
#  define ELEM_SCALAPACK(name) name
# endif
#endif

// TODO: Think of how to better decouple the following components

// Declare the intertwined core parts of our library
#include "elemental/core/Timer/decl.hpp"
#include "elemental/core/Memory.hpp"
#include "elemental/core/Complex/decl.hpp"
#include "elemental/core/types/decl.hpp"
#include "elemental/core/imports/mpi.hpp"
#include "elemental/core/imports/choice.hpp"
#include "elemental/core/imports/mpi_choice.hpp"
#include "elemental/core/environment/decl.hpp"
#include "elemental/core/indexing/decl.hpp"
#include "elemental/core/imports/blas.hpp"
#include "elemental/core/imports/lapack.hpp"
#include "elemental/core/imports/flame.hpp"
#include "elemental/core/imports/pmrrr.hpp"
#include "elemental/core/imports/scalapack.hpp"

#include "elemental/core/Matrix/forward_decl.hpp"
#include "elemental/core/DistMatrix/forward_decl.hpp"
#include "elemental/core/BlockDistMatrix/forward_decl.hpp"
#include "elemental/core/Matrix.hpp"
#include "elemental/core/Grid/decl.hpp"
#include "elemental/core/DistMatrix.hpp"
#include "elemental/core/BlockDistMatrix.hpp"

// Implement the intertwined parts of the library
#include "elemental/core/Timer/impl.hpp"
#include "elemental/core/Complex/impl.hpp"
#include "elemental/core/types/impl.hpp"
#include "elemental/core/Grid/impl.hpp"
#include "elemental/core/environment/impl.hpp"
#include "elemental/core/indexing/impl.hpp"

// Declare and implement the decoupled parts of the core of the library
// (perhaps these should be moved into their own directory?)
#include "elemental/core/views/View.hpp"
#include "elemental/core/views/Partition.hpp"
#include "elemental/core/views/Repartition.hpp"
#include "elemental/core/views/SlidePartition.hpp"
#include "elemental/core/random/decl.hpp"
#include "elemental/core/random/impl.hpp"
#include "elemental/core/AxpyInterface/decl.hpp"
#include "elemental/core/AxpyInterface/impl.hpp"

#endif // ifndef ELEM_CORE_HPP
