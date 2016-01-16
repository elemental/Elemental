/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El.h"
using namespace El;

extern "C" {

#define C_PROTO(SIG,SIGBASE,T) \
  /* View with index ranges */ \
  ElError ElView_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG B, ElRange_i I, ElRange_i J ) \
  { EL_TRY( View( \
      *CReflect(A), *CReflect(B), CReflect(I), CReflect(J) ) ) } \
  ElError ElViewDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B, ElRange_i I, ElRange_i J ) \
  { EL_TRY( View( \
      *CReflect(A), *CReflect(B), CReflect(I), CReflect(J) ) ) } \
  ElError ElLockedView_ ## SIG \
  ( ElMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, ElRange_i I, ElRange_i J ) \
  { EL_TRY( LockedView( \
      *CReflect(A), *CReflect(B), CReflect(I), CReflect(J) ) ) } \
  ElError ElLockedViewDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    ElRange_i I, ElRange_i J ) \
  { EL_TRY( LockedView( \
      *CReflect(A), *CReflect(B), CReflect(I), CReflect(J) ) ) } \
  /* View with an offset and size */ \
  ElError ElViewOffset_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG B, \
    ElInt i, ElInt j, ElInt height, ElInt width ) \
  { EL_TRY( View( *CReflect(A), *CReflect(B), i, j, height, width ) ) } \
  ElError ElViewOffsetDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B, \
    ElInt i, ElInt j, ElInt height, ElInt width ) \
  { EL_TRY( View( *CReflect(A), *CReflect(B), i, j, height, width ) ) } \
  ElError ElLockedViewOffset_ ## SIG \
  ( ElMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    ElInt i, ElInt j, ElInt height, ElInt width ) \
  { EL_TRY( LockedView( *CReflect(A), *CReflect(B), i, j, height, width ) ) } \
  ElError ElLockedViewOffsetDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    ElInt i, ElInt j, ElInt height, ElInt width ) \
  { EL_TRY( LockedView( *CReflect(A), *CReflect(B), i, j, height, width ) ) } \
  /* View a full matrix */ \
  ElError ElViewFull_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( View( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElViewFullDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( View( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElLockedViewFull_ ## SIG \
  ( ElMatrix_ ## SIG A, ElConstMatrix_ ## SIG B ) \
  { EL_TRY( LockedView( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElLockedViewFullDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B ) \
  { EL_TRY( LockedView( *CReflect(A), *CReflect(B) ) ) }

#include "El/macros/CInstantiate.h"

} // extern "C"
