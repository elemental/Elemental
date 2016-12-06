/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include <El.h>
using namespace El;

extern "C" {

#define C_PROTO_BASE(SIG,SIGBASE,T) \
  /* Median
     ====== */ \
  ElError ElMedian_ ## SIG \
  ( ElConstMatrix_ ## SIG x, ElValueInt_ ## SIG *median ) \
  { EL_TRY( *median = CReflect(Median(*CReflect(x))) ) } \
  ElError ElMedianDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG x, ElValueInt_ ## SIG *median ) \
  { EL_TRY( *median = CReflect(Median(*CReflect(x))) ) } \
  /* Sort
     ==== */ \
  ElError ElSort_ ## SIG ( ElMatrix_ ## SIG X, ElSortType sort ) \
  { EL_TRY( Sort( *CReflect(X), CReflect(sort) ) ) } \
  ElError ElSortDist_ ## SIG ( ElDistMatrix_ ## SIG X, ElSortType sort ) \
  { EL_TRY( Sort( *CReflect(X), CReflect(sort) ) ) } \
  ElError ElTaggedSort_ ## SIG \
  ( ElConstMatrix_ ## SIG X, ElSortType sort, \
    ElValueInt_ ## SIG *taggedOrder ) \
  { EL_TRY( \
        auto tagVec = TaggedSort( *CReflect(X), CReflect(sort) ); \
        MemCopy( taggedOrder, CReflect(tagVec.data()), tagVec.size() ); \
    ) } \
  ElError ElTaggedSortDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG X, ElSortType sort, \
    ElValueInt_ ## SIG *taggedOrder ) \
  { EL_TRY( \
        auto tagVec = TaggedSort( *CReflect(X), CReflect(sort) ); \
        MemCopy( taggedOrder, CReflect(tagVec.data()), tagVec.size() ); \
    ) }

#define C_PROTO_INT(SIG,T)     C_PROTO_BASE(SIG,SIG,T)
#define C_PROTO_REAL(SIG,Real) C_PROTO_BASE(SIG,SIG,Real)

#define EL_NO_COMPLEX_PROTO
#include <El/macros/CInstantiate.h>

} // extern "C"
