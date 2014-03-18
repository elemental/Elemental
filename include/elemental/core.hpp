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

#ifdef RELEASE
# define DEBUG_ONLY(cmd) 
#else
# define DEBUG_ONLY(cmd) cmd;
#endif

// If defined, the _OPENMP macro contains the date of the specification
#ifdef HAVE_OPENMP
# include <omp.h>
# define PARALLEL_FOR _Pragma("omp parallel for")
# ifdef HAVE_OMP_COLLAPSE
#  define PARALLEL_FOR_COLLAPSE2 _Pragma("omp parallel for collapse(2)")
# else
#  define PARALLEL_FOR_COLLAPSE2 PARALLEL_FOR
# endif
#else
# define PARALLEL_FOR 
# define PARALLEL_FOR_COLLAPSE2
#endif

#ifdef AVOID_OMP_FMA
# define FMA_PARALLEL_FOR 
#else
# define FMA_PARALLEL_FOR PARALLEL_FOR
#endif
#ifdef PARALLELIZE_INNER_LOOPS
# define INNER_PARALLEL_FOR           PARALLEL_FOR
# define INNER_PARALLEL_FOR_COLLAPSE2 PARALLEL_FOR_COLLAPSE2
# define OUTER_PARALLEL_FOR 
# define OUTER_PARALLEL_FOR_COLLAPSE2
#else
# define INNER_PARALLEL_FOR
# define INNER_PARALLEL_FOR_COLLAPSE2
# define OUTER_PARALLEL_FOR           PARALLEL_FOR
# define OUTER_PARALLEL_FOR_COLLAPSE2 PARALLEL_FOR_COLLAPSE2
#endif

#if defined(BLAS_POST)
#define BLAS(name) name ## _
#else
#define BLAS(name) name
#endif

#if defined(LAPACK_POST)
#define LAPACK(name) name ## _
#else
#define LAPACK(name) name
#endif

// TODO: Think of how to better decouple the following components

// Declare the intertwined core parts of our library
#include "elemental/core/timer/decl.hpp"
#include "elemental/core/memory/decl.hpp"
#include "elemental/core/complex/decl.hpp"
#include "elemental/core/types/decl.hpp"
#include "elemental/core/matrix/forward_decl.hpp"
#include "elemental/core/dist_matrix/forward_decl.hpp"
//#include "elemental/core/block_dist_matrix/forward_decl.hpp"
#include "elemental/core/view/decl.hpp"
#include "elemental/core/matrix.hpp"
#include "elemental/core/imports/mpi.hpp"
#include "elemental/core/grid/decl.hpp"
#include "elemental/core/dist_matrix.hpp"
//#include "elemental/core/block_dist_matrix.hpp"
#include "elemental/core/imports/choice.hpp"
#include "elemental/core/imports/mpi_choice.hpp"
#include "elemental/core/environment/decl.hpp"
#include "elemental/core/indexing/decl.hpp"
#include "elemental/core/imports/blas.hpp"
#include "elemental/core/imports/lapack.hpp"
#include "elemental/core/imports/flame.hpp"
#include "elemental/core/imports/pmrrr.hpp"

// Implement the intertwined parts of the library
#include "elemental/core/timer/impl.hpp"
#include "elemental/core/memory/impl.hpp"
#include "elemental/core/complex/impl.hpp"
#include "elemental/core/types/impl.hpp"
#include "elemental/core/grid/impl.hpp"
#include "elemental/core/environment/impl.hpp"
#include "elemental/core/indexing/impl.hpp"

// Declare and implement the decoupled parts of the core of the library
// (perhaps these should be moved into their own directory?)
#include "elemental/core/view/impl.hpp"
#include "elemental/core/partition/decl.hpp"
#include "elemental/core/partition/impl.hpp"
#include "elemental/core/repartition/decl.hpp"
#include "elemental/core/repartition/impl.hpp"
#include "elemental/core/slide_partition/decl.hpp"
#include "elemental/core/slide_partition/impl.hpp"
#include "elemental/core/random/decl.hpp"
#include "elemental/core/random/impl.hpp"
#include "elemental/core/axpy_interface/decl.hpp"
#include "elemental/core/axpy_interface/impl.hpp"

#endif // ifndef ELEM_CORE_HPP
