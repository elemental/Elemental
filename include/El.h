/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_C_H
#define EL_C_H

#include "mpi.h"

/* TODO: A better include structure for the C interface */
#include "El/config.h"
#ifdef EL_HAVE_F90_INTERFACE
# include "El/FCMangle.h"
#endif

#include "El/core/types.h"
#include "El/core/environment.h"
#include "El/core/imports/mpi.h"
#include "El/core/Grid.h"
#include "El/core/Matrix.h"
#include "El/core/DistMatrix.h"
#include "El/core/Graph.h"
#include "El/core/DistGraph.h"
#include "El/core/SparseMatrix.h"
#include "El/core/DistSparseMatrix.h"
#include "El/core/DistMultiVec.h"
#include "El/core/View.h"
#include "El/core/FlamePart/Merge.h"
#include "El/core/FlamePart/Partition.h"
#include "El/core/FlamePart/Repartition.h"
#include "El/core/FlamePart/SlidePartition.h"

#include "El/io.h"

#include "El/blas_like/level1.h"
#include "El/blas_like/level2.h"
#include "El/blas_like/level3.h"

#include "El/lapack_like/perm.h"
#include "El/lapack_like/reflect.h"
#include "El/lapack_like/util.h"

#include "El/lapack_like/condense.h"
#include "El/lapack_like/factor.h"

#include "El/lapack_like/spectral.h"
#include "El/lapack_like/funcs.h"

#include "El/lapack_like/solve.h"
#include "El/lapack_like/euclidean_min.h"

#include "El/lapack_like/props.h"

#include "El/matrices.h"

#include "El/control.h"
#include "El/optimization/solvers.h"
#include "El/optimization/models.h"
#include "El/optimization/prox.h"
#include "El/optimization/util.h"

#include "El/lattice.h"

#ifdef __cplusplus
#include "El/CReflect.hpp"
#endif

#endif /* ifndef EL_C_H */
