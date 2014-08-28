/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_UTIL_C_H
#define EL_UTIL_C_H

#include "El/core/DistMatrix-C.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Median
   ====== */
ElError ElMedian_i( ElConstMatrix_i x, ElValueInt_i* median );
ElError ElMedian_s( ElConstMatrix_s x, ElValueInt_s* median );
ElError ElMedian_d( ElConstMatrix_d x, ElValueInt_d* median );

ElError ElMedianDist_i( ElConstDistMatrix_i x, ElValueInt_i* median );
ElError ElMedianDist_s( ElConstDistMatrix_s x, ElValueInt_s* median );
ElError ElMedianDist_d( ElConstDistMatrix_d x, ElValueInt_d* median );

/* Sort
   ==== */
ElError ElSort_i( ElMatrix_i X, ElSortType sort );
ElError ElSort_s( ElMatrix_s X, ElSortType sort );
ElError ElSort_d( ElMatrix_d X, ElSortType sort );

ElError ElSortDist_i( ElDistMatrix_i X, ElSortType sort );
ElError ElSortDist_s( ElDistMatrix_s X, ElSortType sort );
ElError ElSortDist_d( ElDistMatrix_d X, ElSortType sort );

ElError ElTaggedSort_i
( ElConstMatrix_i x, ElSortType sort, ElValueInt_i* taggedOrder );
ElError ElTaggedSort_s
( ElConstMatrix_s x, ElSortType sort, ElValueInt_s* taggedOrder );
ElError ElTaggedSort_d
( ElConstMatrix_d x, ElSortType sort, ElValueInt_d* taggedOrder );

ElError ElTaggedSortDist_i
( ElConstDistMatrix_i x, ElSortType sort, ElValueInt_i* taggedOrder );
ElError ElTaggedSortDist_s
( ElConstDistMatrix_s x, ElSortType sort, ElValueInt_s* taggedOrder );
ElError ElTaggedSortDist_d
( ElConstDistMatrix_d x, ElSortType sort, ElValueInt_d* taggedOrder );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_UTIL_C_H */
