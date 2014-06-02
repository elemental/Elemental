/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CORE_HPP
#define EL_CORE_HPP

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
#ifdef EL_RELEASE
# define DEBUG_ONLY(cmd) 
#else
# define DEBUG_ONLY(cmd) cmd;
#endif

// If defined, the _OPENMP macro contains the date of the specification
#ifdef EL_HAVE_OPENMP
# include <omp.h>
# define EL_PARALLEL_FOR _Pragma("omp parallel for")
# ifdef EL_HAVE_OMP_COLLAPSE
#  define EL_PARALLEL_FOR_COLLAPSE2 _Pragma("omp parallel for collapse(2)")
# else
#  define EL_PARALLEL_FOR_COLLAPSE2 EL_PARALLEL_FOR
# endif
#else
# define EL_PARALLEL_FOR 
# define EL_PARALLEL_FOR_COLLAPSE2
#endif

#ifdef EL_AVOID_OMP_FMA
# define EL_FMA_PARALLEL_FOR 
#else
# define EL_FMA_PARALLEL_FOR EL_PARALLEL_FOR
#endif
#ifdef EL_PARALLELIZE_INNER_LOOPS
# define EL_INNER_PARALLEL_FOR           EL_PARALLEL_FOR
# define EL_INNER_PARALLEL_FOR_COLLAPSE2 EL_PARALLEL_FOR_COLLAPSE2
# define EL_OUTER_PARALLEL_FOR 
# define EL_OUTER_PARALLEL_FOR_COLLAPSE2
#else
# define EL_INNER_PARALLEL_FOR
# define EL_INNER_PARALLEL_FOR_COLLAPSE2
# define EL_OUTER_PARALLEL_FOR           EL_PARALLEL_FOR
# define EL_OUTER_PARALLEL_FOR_COLLAPSE2 EL_PARALLEL_FOR_COLLAPSE2
#endif

#ifdef EL_HAVE_NOEXCEPT
# define EL_NOEXCEPT noexcept
#else
# define EL_NOEXCEPT
#endif

#if defined(EL_BLAS_POST)
# define EL_BLAS(name) name ## _
#else
# define EL_BLAS(name) name
#endif

#if defined(EL_LAPACK_POST)
# define EL_LAPACK(name) name ## _
#else
# define EL_LAPACK(name) name
#endif

#if defined(EL_HAVE_SCALAPACK)
# if defined(EL_SCALAPACK_POST)
#  define EL_SCALAPACK(name) name ## _
# else
#  define EL_SCALAPACK(name) name
# endif
#endif

// TODO: Think of how to better decouple the following components

// Declare the intertwined core parts of our library
#include "El/core/Timer.hpp"
#include "El/core/Memory.hpp"
#include "El/core/Scalar/decl.hpp"
#include "El/core/types.hpp"
#include "El/core/imports/mpi.hpp"
#include "El/core/imports/choice.hpp"
#include "El/core/imports/mpi_choice.hpp"
#include "El/core/environment/decl.hpp"
#include "El/core/indexing/decl.hpp"
#include "El/core/imports/blas.hpp"
#include "El/core/imports/lapack.hpp"
#include "El/core/imports/flame.hpp"
#include "El/core/imports/pmrrr.hpp"
#include "El/core/imports/scalapack.hpp"

namespace El {

template<typename T> class Matrix;

template<typename T> class AbstractDistMatrix;
template<typename T> class AbstractBlockDistMatrix;

template<typename T,Dist U=MC,Dist V=MR> class GeneralDistMatrix;
template<typename T,Dist U=MC,Dist V=MR> class GeneralBlockDistMatrix;

template<typename T,Dist U=MC,Dist V=MR> class DistMatrix;
template<typename T,Dist U=MC,Dist V=MR> class BlockDistMatrix;

} // namespace El

#include "El/core/Matrix.hpp"
#include "El/core/Grid.hpp"
#include "El/core/DistMatrix.hpp"
#include "El/core/BlockDistMatrix.hpp"

// Implement the intertwined parts of the library
#include "El/core/Scalar/impl.hpp"
#include "El/core/environment/impl.hpp"
#include "El/core/indexing/impl.hpp"

// Declare and implement the decoupled parts of the core of the library
// (perhaps these should be moved into their own directory?)
#include "El/core/views/View.hpp"
#include "El/core/views/Partition.hpp"
#include "El/core/views/Repartition.hpp"
#include "El/core/views/SlidePartition.hpp"
#include "El/core/random/decl.hpp"
#include "El/core/random/impl.hpp"
#include "El/core/AxpyInterface.hpp"

#endif // ifndef EL_CORE_HPP
