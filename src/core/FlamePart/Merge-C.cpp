/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El/core/FlamePart.hpp"

#include "El.h"
using namespace El;

extern "C" {

#define C_PROTO(SIG,SIGBASE,T) \
  /* Horizontally merge two contiguous matrices */ \
  ElError ElMerge1x2_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG BL, ElMatrix_ ## SIG BR ) \
  { EL_TRY( Merge1x2( *CReflect(A), *CReflect(BL), *CReflect(BR) ) ) } \
  ElError ElMerge1x2Dist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG BL, ElDistMatrix_ ## SIG BR ) \
  { EL_TRY( Merge1x2( *CReflect(A), *CReflect(BL), *CReflect(BR) ) ) } \
  ElError ElLockedMerge1x2_ ## SIG \
  ( ElMatrix_ ## SIG A, ElConstMatrix_ ## SIG BL, ElConstMatrix_ ## SIG BR ) \
  { EL_TRY( LockedMerge1x2( *CReflect(A), *CReflect(BL), *CReflect(BR) ) ) } \
  ElError ElLockedMerge1x2Dist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG BL, ElConstDistMatrix_ ## SIG BR ) \
  { EL_TRY( LockedMerge1x2( *CReflect(A), *CReflect(BL), *CReflect(BR) ) ) } \
  /* Vertically merge two contiguous matrices */ \
  ElError ElMerge2x1_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG BL, ElMatrix_ ## SIG BR ) \
  { EL_TRY( Merge2x1( *CReflect(A), *CReflect(BL), *CReflect(BR) ) ) } \
  ElError ElMerge2x1Dist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG BL, ElDistMatrix_ ## SIG BR ) \
  { EL_TRY( Merge2x1( *CReflect(A), *CReflect(BL), *CReflect(BR) ) ) } \
  ElError ElLockedMerge2x1_ ## SIG \
  ( ElMatrix_ ## SIG A, ElConstMatrix_ ## SIG BL, ElConstMatrix_ ## SIG BR ) \
  { EL_TRY( LockedMerge2x1( *CReflect(A), *CReflect(BL), *CReflect(BR) ) ) } \
  ElError ElLockedMerge2x1Dist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG BL, ElConstDistMatrix_ ## SIG BR ) \
  { EL_TRY( LockedMerge2x1( *CReflect(A), *CReflect(BL), *CReflect(BR) ) ) } \
  /* Merge a contiguous 2x2 matrix of blocks */ \
  ElError ElMerge2x2_ ## SIG \
  ( ElMatrix_ ## SIG A, \
    ElMatrix_ ## SIG BTL, ElMatrix_ ## SIG BTR, \
    ElMatrix_ ## SIG BBL, ElMatrix_ ## SIG BBR ) \
  { EL_TRY( Merge2x2( *CReflect(A), \
      *CReflect(BTL), *CReflect(BTR), *CReflect(BBL), *CReflect(BBR)) ) } \
  ElError ElMerge2x2Dist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG BTL, ElDistMatrix_ ## SIG BTR, \
    ElDistMatrix_ ## SIG BBL, ElDistMatrix_ ## SIG BBR ) \
  { EL_TRY( Merge2x2( *CReflect(A), \
      *CReflect(BTL), *CReflect(BTR), *CReflect(BBL), *CReflect(BBR)) ) } \
  ElError ElLockedMerge2x2_ ## SIG \
  ( ElMatrix_ ## SIG A, \
    ElConstMatrix_ ## SIG BTL, ElConstMatrix_ ## SIG BTR, \
    ElConstMatrix_ ## SIG BBL, ElConstMatrix_ ## SIG BBR ) \
  { EL_TRY( LockedMerge2x2( *CReflect(A), \
      *CReflect(BTL), *CReflect(BTR), *CReflect(BBL), *CReflect(BBR)) ) } \
  ElError ElLockedMerge2x2Dist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElConstDistMatrix_ ## SIG BTL, ElConstDistMatrix_ ## SIG BTR, \
    ElConstDistMatrix_ ## SIG BBL, ElConstDistMatrix_ ## SIG BBR ) \
  { EL_TRY( LockedMerge2x2( *CReflect(A), \
      *CReflect(BTL), *CReflect(BTR), *CReflect(BBL), *CReflect(BBR)) ) }

#include "El/macros/CInstantiate.h"

} // extern "C"
