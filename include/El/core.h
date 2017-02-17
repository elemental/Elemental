/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CORE_C_H
#define EL_CORE_C_H

#include <mpi.h>

/* TODO: A better include structure for the C interface */
#include <El/config.h>
#ifdef EL_HAVE_F90_INTERFACE
# include <El/FCMangle.h>
#endif

#include <El/core/types.h>
#include <El/core/environment.h>
#include <El/core/imports/mpi.h>
#include <El/core/Grid.h>
#include <El/core/Matrix.h>
#include <El/core/DistMatrix.h>
#include <El/core/Graph.h>
#include <El/core/DistGraph.h>
#include <El/core/SparseMatrix.h>
#include <El/core/DistSparseMatrix.h>
#include <El/core/DistMultiVec.h>
#include <El/core/Permutation.h>

#include <El/core/View.h>
#include <El/core/FlamePart/Merge.h>
#include <El/core/FlamePart/Partition.h>
#include <El/core/FlamePart/Repartition.h>
#include <El/core/FlamePart/SlidePartition.h>

#ifdef __cplusplus
#include <El/core/CReflect.hpp>
#endif

#endif /* ifndef EL_CORE_C_H */
