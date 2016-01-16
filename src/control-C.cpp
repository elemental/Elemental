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
  /* Lyapunov 
     ======== */ \
  ElError ElLyapunov_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG C, \
    ElMatrix_ ## SIG X ) \
  { EL_TRY( Lyapunov( *CReflect(A), *CReflect(C), *CReflect(X) ) ) } \
  ElError ElLyapunovDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG C, \
    ElDistMatrix_ ## SIG X ) \
  { EL_TRY( Lyapunov( *CReflect(A), *CReflect(C), *CReflect(X) ) ) } \
  /* Ricatti
     ======= */ \
  ElError ElRicatti_ ## SIG \
  ( ElUpperOrLower uplo, \
    ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG K, \
    ElConstMatrix_ ## SIG L, ElMatrix_ ## SIG X ) \
  { EL_TRY( Ricatti( CReflect(uplo), \
      *CReflect(A), *CReflect(K), *CReflect(L), *CReflect(X) ) ) } \
  ElError ElRicattiDist_ ## SIG \
  ( ElUpperOrLower uplo, \
    ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG K, \
    ElConstDistMatrix_ ## SIG L, ElDistMatrix_ ## SIG X ) \
  { EL_TRY( Ricatti( CReflect(uplo), \
      *CReflect(A), *CReflect(K), *CReflect(L), *CReflect(X) ) ) } \
  /* Preformed
     --------- */ \
  ElError ElRicattiPreformed_ ## SIG \
  ( ElMatrix_ ## SIG W, ElMatrix_ ## SIG X ) \
  { EL_TRY( Ricatti( *CReflect(W), *CReflect(X) ) ) } \
  ElError ElRicattiPreformedDist_ ## SIG \
  ( ElDistMatrix_ ## SIG W, ElDistMatrix_ ## SIG X ) \
  { EL_TRY( Ricatti( *CReflect(W), *CReflect(X) ) ) } \
  /* Sylvester
     ========= */ \
  ElError ElSylvester_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    ElConstMatrix_ ## SIG C, ElMatrix_ ## SIG X ) \
  { EL_TRY( Sylvester( \
      *CReflect(A), *CReflect(B), *CReflect(C), *CReflect(X) ) ) } \
  ElError ElSylvesterDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    ElConstDistMatrix_ ## SIG C, ElDistMatrix_ ## SIG X ) \
  { EL_TRY( Sylvester( \
      *CReflect(A), *CReflect(B), *CReflect(C), *CReflect(X) ) ) } \
  /* Preformed
     --------- */ \
  ElError ElSylvesterPreformed_ ## SIG \
  ( ElInt m, ElMatrix_ ## SIG W, ElMatrix_ ## SIG X ) \
  { EL_TRY( Sylvester( m, *CReflect(W), *CReflect(X) ) ) } \
  ElError ElSylvesterPreformedDist_ ## SIG \
  ( ElInt m, ElDistMatrix_ ## SIG W, ElDistMatrix_ ## SIG X ) \
  { EL_TRY( Sylvester( m, *CReflect(W), *CReflect(X) ) ) }

#define EL_NO_INT_PROTO
#include "El/macros/CInstantiate.h"

} // extern "C"
