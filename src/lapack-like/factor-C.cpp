/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El-C.h"
using namespace El;

#define DM_CAST(T,A) dynamic_cast<DistMatrix<T>&>(*Reinterpret(A))
#define DM_CAST_CONST(T,A) dynamic_cast<const DistMatrix<T>&>(*Reinterpret(A))

#define DM_VC_STAR_CAST(T,A) \
  dynamic_cast<DistMatrix<T,VC,STAR>&>(*Reinterpret(A))
#define DM_VC_STAR_CAST_CONST(T,A) \
  dynamic_cast<const DistMatrix<T,VC,STAR>&>(*Reinterpret(A))

extern "C" {

#define C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* Cholesky (no pivoting) */ \
  ElError ElCholesky_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( Cholesky( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  ElError ElCholeskyDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Cholesky( Reinterpret(uplo), DM_CAST(F,A) ) ) } \
  /* Reverse Cholesky (no pivoting) */ \
  ElError ElReverseCholesky_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( ReverseCholesky( Reinterpret(uplo), *Reinterpret(A) ) ) } \
  ElError ElReverseCholeskyDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( ReverseCholesky( Reinterpret(uplo), DM_CAST(F,A) ) ) } \
  /* Cholesky (full pivoting) */ \
  ElError ElCholeskyPiv_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElMatrix_i p ) \
  { EL_TRY( \
      Cholesky( Reinterpret(uplo), *Reinterpret(A), *Reinterpret(p) ) ) } \
  ElError ElCholeskyPivDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElDistMatrix_i p ) \
  { EL_TRY( \
      Cholesky( Reinterpret(uplo), DM_CAST(F,A), DM_VC_STAR_CAST(Int,p) ) ) } \
  /* Cholesky low-rank modification */ \
  ElError ElCholeskyMod_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG T, \
    Base<F> alpha, ElMatrix_ ## SIG V ) \
  { EL_TRY( \
      CholeskyMod( \
        Reinterpret(uplo), *Reinterpret(T), alpha, *Reinterpret(V) ) ) } \
  ElError ElCholeskyModDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG T, \
    Base<F> alpha, ElDistMatrix_ ## SIG V ) \
  { EL_TRY( \
      CholeskyMod( \
        Reinterpret(uplo), DM_CAST(F,T), alpha, DM_CAST(F,V) ) ) }

#define C_PROTO(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F)

#define C_PROTO_COMPLEX(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F)

#define EL_NO_INT_PROTO
#include "El/macros/CInstantiate.h"

} // extern "C"
