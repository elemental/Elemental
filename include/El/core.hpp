/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CORE_HPP
#define EL_CORE_HPP

// This would ideally be included within core/imports/mpi.hpp, but it is 
// well-known that this must often be included first.
#include <mpi.h>

#include <algorithm>
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
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <random>
#include <type_traits>
#include <vector>
#include <iomanip>

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

// Declare the intertwined core parts of our library
#include "El/core/imports/valgrind.hpp"
#include "El/core/imports/omp.hpp"
#include "El/core/imports/mpc.hpp"
#include "El/core/Memory.hpp"
#include "El/core/Element/decl.hpp"
#include "El/core/types.hpp"
#include "El/core/imports/mpi.hpp"
#include "El/core/imports/choice.hpp"
#include "El/core/imports/mpi_choice.hpp"
#include "El/core/environment/decl.hpp"

#include "El/core/Timer.hpp"
#include "El/core/indexing/decl.hpp"
#include "El/core/imports/blas.hpp"
#include "El/core/imports/lapack.hpp"
#include "El/core/imports/flame.hpp"
#include "El/core/imports/mkl.hpp"
#include "El/core/imports/openblas.hpp"
#include "El/core/imports/pmrrr.hpp"
#include "El/core/imports/scalapack.hpp"

namespace El {

template<typename T=double> class Matrix;

template<typename T=double> class AbstractDistMatrix;

template<typename T=double> class ElementalMatrix;
template<typename T=double> class BlockMatrix;

template<typename T=double,Dist U=MC,Dist V=MR,DistWrap wrap=ELEMENT>
class DistMatrix;

} // namespace El

#include "El/core/Matrix.hpp"
#include "El/core/Grid.hpp"
#include "El/core/DistMatrix.hpp"
#include "El/core/Proxy.hpp"

// Implement the intertwined parts of the library
#include "El/core/Element/impl.hpp"
#include "El/core/environment/impl.hpp"
#include "El/core/indexing/impl.hpp"

// Declare and implement the decoupled parts of the core of the library
// (perhaps these should be moved into their own directory?)
#include "El/core/View.hpp"
#include "El/core/flame_part/Merge.hpp"
#include "El/core/flame_part/Partition.hpp"
#include "El/core/flame_part/Repartition.hpp"
#include "El/core/flame_part/SlidePartition.hpp"
#include "El/core/random/decl.hpp"
#include "El/core/random/impl.hpp"

#include "El/core/Graph.hpp"
// TODO: Sequential map
//#include "El/core/Map.hpp"
#include "El/core/SparseMatrix.hpp"

#include "El/core/DistGraph.hpp"
#include "El/core/DistMap.hpp"
#include "El/core/DistMultiVec.hpp"
#include "El/core/DistSparseMatrix.hpp"

#endif // ifndef EL_CORE_HPP
