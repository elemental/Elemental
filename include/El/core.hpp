/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CORE_HPP
#define EL_CORE_HPP

// This would ideally be included within core/imports/mpi.hpp, but it is 
// well-known that this must often be included first.
#include <mpi.h>

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <random>
#include <type_traits> // std::enable_if
#include <vector>

// The DEBUG_ONLY and RELEASE_ONLY macros are, to the best of my knowledge, 
// the only preprocessor names defined by Elemental that is not namespaced 
// with "EL". Given how frequently they are used, I will leave it as-is 
// unless/until a user/developer complains.
#ifdef EL_RELEASE
# define DEBUG_ONLY(cmd) 
# define RELEASE_ONLY(cmd) cmd;
#else
# define DEBUG_ONLY(cmd) cmd;
# define RELEASE_ONLY(cmd)
#endif

#ifdef EL_HAVE_NO_EXCEPT
# define EL_NO_EXCEPT noexcept
#else
# define EL_NO_EXCEPT
#endif

#ifdef EL_RELEASE
# define EL_NO_RELEASE_EXCEPT EL_NO_EXCEPT
#else
# define EL_NO_RELEASE_EXCEPT
#endif

#define EL_CONCAT2(name1,name2) name1 ## name2
#define EL_CONCAT(name1,name2) EL_CONCAT2(name1,name2)

// TODO: Think of how to better decouple the following components

#ifdef EL_HAVE_QUADMATH
#include <quadmath.h>
#endif

namespace El {

typedef unsigned char byte;

// If these are changes, you must make sure that they have 
// existing MPI datatypes. This is only sometimes true for 'long long'
#ifdef EL_USE_64BIT_INTS
typedef long long int Int;
typedef long long unsigned Unsigned;
#else
typedef int Int;
typedef unsigned Unsigned;
#endif

#ifdef EL_HAVE_QUAD
typedef __float128 Quad;
#endif

} // namespace El

// Declare the intertwined core parts of our library
#include <El/core/imports/valgrind.hpp>
#include <El/core/imports/omp.hpp>
#include <El/core/imports/qd.hpp>
#include <El/core/imports/mpc.hpp>

#include <El/core/Element/decl.hpp>
#include <El/core/types.hpp>
#include <El/core/Serialize.hpp>
#include <El/core/imports/mpi.hpp>
#include <El/core/imports/choice.hpp>
#include <El/core/imports/mpi_choice.hpp>
#include <El/core/environment/decl.hpp>

#include <El/core/Timer.hpp>
#include <El/core/indexing/decl.hpp>
#include <El/core/imports/blas.hpp>
#include <El/core/imports/lapack.hpp>
#include <El/core/imports/flame.hpp>
#include <El/core/imports/mkl.hpp>
#include <El/core/imports/openblas.hpp>
#include <El/core/imports/pmrrr.hpp>
#include <El/core/imports/scalapack.hpp>

#include <El/core/limits.hpp>

#include <El/core/Memory.hpp>

namespace El {

template<typename T=double> class Matrix;

template<typename T=double> class AbstractDistMatrix;

template<typename T=double> class ElementalMatrix;
template<typename T=double> class BlockMatrix;

template<typename T=double,Dist U=MC,Dist V=MR,DistWrap wrap=ELEMENT>
class DistMatrix;

} // namespace El

#include <El/core/Matrix/decl.hpp>
#include <El/core/Graph/decl.hpp>
#include <El/core/DistMap/decl.hpp>
#include <El/core/DistGraph/decl.hpp>
#include <El/core/SparseMatrix/decl.hpp>
#include <El/core/DistSparseMatrix/decl.hpp>
#include <El/core/DistMultiVec/decl.hpp>
#include <El/core/View/decl.hpp>
#include <El/blas_like/level1/decl.hpp>

#include <El/core/Matrix/impl.hpp>
#include <El/core/Grid.hpp>
#include <El/core/DistMatrix.hpp>
#include <El/core/Proxy.hpp>

// Implement the intertwined parts of the library
#include <El/core/Element/impl.hpp>
#include <El/core/environment/impl.hpp>
#include <El/core/indexing/impl.hpp>

// Declare and implement the decoupled parts of the core of the library
// (perhaps these should be moved into their own directory?)
#include <El/core/View/impl.hpp>
#include <El/core/FlamePart.hpp>
#include <El/core/random/decl.hpp>
#include <El/core/random/impl.hpp>

// TODO: Sequential map
//#include <El/core/Map.hpp>
#include <El/core/SparseMatrix/impl.hpp>

#include <El/core/DistMap.hpp>
#include <El/core/DistMultiVec/impl.hpp>
#include <El/core/DistSparseMatrix/impl.hpp>

#endif // ifndef EL_CORE_HPP
