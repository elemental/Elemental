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

ElError ElSymvCtrlDefault_s( ElSymvCtrl* ctrl )
{
    ctrl->bsize = LocalSymvBlocksize<float>();
    ctrl->avoidTrmvBasedLocalSymv = true;
    return EL_SUCCESS;
}
ElError ElSymvCtrlDefault_d( ElSymvCtrl* ctrl )
{
    ctrl->bsize = LocalSymvBlocksize<double>();
    ctrl->avoidTrmvBasedLocalSymv = true;
    return EL_SUCCESS;
}
ElError ElSymvCtrlDefault_c( ElSymvCtrl* ctrl )
{
    ctrl->bsize = LocalSymvBlocksize<Complex<float>>();
    ctrl->avoidTrmvBasedLocalSymv = true;
    return EL_SUCCESS;
}
ElError ElSymvCtrlDefault_z( ElSymvCtrl* ctrl )
{
    ctrl->bsize = LocalSymvBlocksize<Complex<double>>();
    ctrl->avoidTrmvBasedLocalSymv = true;
    return EL_SUCCESS;
}

#define C_PROTO_BASE(SIG,SIGBASE,T) \
  /* Gemv */ \
  ElError ElGemv_ ## SIG \
  ( ElOrientation orientation, \
    CREFLECT(T) alpha, ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG x, \
    CREFLECT(T) beta, ElMatrix_ ## SIG y ) \
  { EL_TRY( \
      Gemv( CReflect(orientation), \
            CReflect(alpha), *CReflect(A), *CReflect(x), \
            CReflect(beta), *CReflect(y) ) ) } \
  ElError ElGemvDist_ ## SIG \
  ( ElOrientation orientation, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG A, \
                       ElConstDistMatrix_ ## SIG x, \
    CREFLECT(T) beta,  ElDistMatrix_ ## SIG y ) \
  { EL_TRY( \
      Gemv( CReflect(orientation), \
            CReflect(alpha), *CReflect(A), *CReflect(x), \
            CReflect(beta), *CReflect(y) ) ) }

#define C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* Ger */ \
  ElError ElGer_ ## SIG \
  ( CREFLECT(F) alpha, ElConstMatrix_ ## SIG x, ElConstMatrix_ ## SIG y, \
    ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      Ger( CReflect(alpha), *CReflect(x), *CReflect(y), \
           *CReflect(A) ) ) } \
  ElError ElGerDist_ ## SIG \
  ( CREFLECT(F) alpha, ElConstDistMatrix_ ## SIG x, \
                       ElConstDistMatrix_ ## SIG y, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      Ger( CReflect(alpha), *CReflect(x), *CReflect(y), *CReflect(A) ) ) } \
  /* Symv */ \
  ElError ElSymv_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(F) alpha, ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG x, \
    CREFLECT(F) beta, ElMatrix_ ## SIG y ) \
  { EL_TRY( \
      Symv( CReflect(uplo), \
            CReflect(alpha), *CReflect(A), *CReflect(x), \
            CReflect(beta), *CReflect(y) ) ) } \
  ElError ElSymvDist_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(F) alpha, ElConstDistMatrix_ ## SIG A, \
                       ElConstDistMatrix_ ## SIG x, \
    CREFLECT(F) beta,  ElDistMatrix_ ## SIG y ) \
  { EL_TRY( \
      Symv( CReflect(uplo), \
            CReflect(alpha), *CReflect(A), *CReflect(x), \
            CReflect(beta), *CReflect(y) ) ) } \
  /* Syr */ \
  ElError ElSyr_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(F) alpha, ElConstMatrix_ ## SIG x, ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      Syr( CReflect(uplo), \
            CReflect(alpha), *CReflect(x), *CReflect(A) ) ) } \
  ElError ElSyrDist_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(F) alpha, ElConstDistMatrix_ ## SIG x, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      Syr( CReflect(uplo), \
           CReflect(alpha), *CReflect(x), *CReflect(A) ) ) } \
  /* Syr2 */ \
  ElError ElSyr2_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(F) alpha, ElConstMatrix_ ## SIG x, ElConstMatrix_ ## SIG y, \
    ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      Syr2( CReflect(uplo), \
            CReflect(alpha), *CReflect(x), *CReflect(y), \
            *CReflect(A) ) ) } \
  ElError ElSyr2Dist_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(F) alpha, ElConstDistMatrix_ ## SIG x, \
                       ElConstDistMatrix_ ## SIG y, \
                       ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      Syr2( CReflect(uplo), \
            CReflect(alpha), *CReflect(x), *CReflect(y), *CReflect(A) ) ) } \
  /* QuasiTrsv */ \
  ElError ElQuasiTrsv_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, ElConstMatrix_ ## SIG A, \
    ElMatrix_ ## SIG x ) \
  { EL_TRY( \
      QuasiTrsv( \
        CReflect(uplo), CReflect(orientation), \
        *CReflect(A), *CReflect(x) ) ) } \
  ElError ElQuasiTrsvDist_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG x ) \
  { EL_TRY( \
      QuasiTrsv( \
        CReflect(uplo), CReflect(orientation), \
        *CReflect(A), *CReflect(x) ) ) } \
  /* Trmv */ \
  ElError ElTrmv_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, ElUnitOrNonUnit diag, \
    ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG x ) \
  { EL_TRY( \
      Trmv( \
        CReflect(uplo), CReflect(orientation), CReflect(diag), \
        *CReflect(A), *CReflect(x) ) ) } \
  /* NOTE: Distributed Trmv not implemented, use Trmm */ \
  /* Trsv */ \
  ElError ElTrsv_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, ElUnitOrNonUnit diag, \
    ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG x ) \
  { EL_TRY( \
      Trsv( \
        CReflect(uplo), CReflect(orientation), CReflect(diag), \
        *CReflect(A), *CReflect(x) ) ) } \
  ElError ElTrsvDist_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, ElUnitOrNonUnit diag, \
    ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG x ) \
  { EL_TRY( \
      Trsv( \
        CReflect(uplo), CReflect(orientation), CReflect(diag), \
        *CReflect(A), *CReflect(x) ) ) }

#define C_PROTO_INT(SIG,T) \
  C_PROTO_BASE(SIG,SIG,T) \
  /* Trr */ \
  ElError ElTrr_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstMatrix_ ## SIG x, ElConstMatrix_ ## SIG y, \
    ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      Trr( CReflect(uplo), \
           CReflect(alpha), *CReflect(x), *CReflect(y), \
           *CReflect(A) ) ) } \
  ElError ElTrrDist_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG x, \
                       ElConstDistMatrix_ ## SIG y, \
                       ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      Trr( CReflect(uplo), \
           CReflect(alpha), *CReflect(x), *CReflect(y), *CReflect(A) ) ) } \
  /* Trr2 */ \
  ElError ElTrr2_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstMatrix_ ## SIG X, ElConstMatrix_ ## SIG Y, \
    ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      Trr2( CReflect(uplo), \
            CReflect(alpha), *CReflect(X), *CReflect(Y), *CReflect(A) ) ) } \
  ElError ElTrr2Dist_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG X, \
                       ElConstDistMatrix_ ## SIG Y, \
                       ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      Trr2( CReflect(uplo), \
            CReflect(alpha), *CReflect(X), *CReflect(Y), *CReflect(A) ) ) } 

#define C_PROTO_REAL(SIG,Real) \
  C_PROTO_BASE(SIG,SIG,Real) \
  C_PROTO_FIELD(SIG,SIG,Real) \
  /* Trr */ \
  ElError ElTrr_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(Real) alpha, ElConstMatrix_ ## SIG x, ElConstMatrix_ ## SIG y, \
    ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      Trr( CReflect(uplo), \
           CReflect(alpha), *CReflect(x), *CReflect(y), \
           *CReflect(A) ) ) } \
  ElError ElTrrDist_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(Real) alpha, ElConstDistMatrix_ ## SIG x, \
                          ElConstDistMatrix_ ## SIG y, \
                          ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      Trr( CReflect(uplo), \
           CReflect(alpha), *CReflect(x), *CReflect(y), *CReflect(A) ) ) } \
  /* Trr2 */ \
  ElError ElTrr2_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(Real) alpha, ElConstMatrix_ ## SIG X, ElConstMatrix_ ## SIG Y, \
    ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      Trr2( CReflect(uplo), \
            CReflect(alpha), *CReflect(X), *CReflect(Y), *CReflect(A) ) ) } \
  ElError ElTrr2Dist_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(Real) alpha, ElConstDistMatrix_ ## SIG X, \
                          ElConstDistMatrix_ ## SIG Y, \
                          ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      Trr2( CReflect(uplo), \
            CReflect(alpha), *CReflect(X), *CReflect(Y), *CReflect(A) ) ) } 

#define C_PROTO_COMPLEX(SIG,SIGBASE,F) \
  C_PROTO_BASE(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* Geru */ \
  ElError ElGeru_ ## SIG \
  ( CREFLECT(F) alpha, ElConstMatrix_ ## SIG x, ElConstMatrix_ ## SIG y, \
    ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      Geru( CReflect(alpha), *CReflect(x), *CReflect(y), \
           *CReflect(A) ) ) } \
  ElError ElGeruDist_ ## SIG \
  ( CREFLECT(F) alpha, ElConstDistMatrix_ ## SIG x, \
                       ElConstDistMatrix_ ## SIG y, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      Geru( CReflect(alpha), *CReflect(x), *CReflect(y), *CReflect(A) ) ) } \
  /* Hemv */ \
  ElError ElHemv_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(F) alpha, ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG x, \
    CREFLECT(F) beta, ElMatrix_ ## SIG y ) \
  { EL_TRY( \
      Hemv( CReflect(uplo), \
            CReflect(alpha), *CReflect(A), *CReflect(x), \
            CReflect(beta), *CReflect(y) ) ) } \
  ElError ElHemvDist_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(F) alpha, ElConstDistMatrix_ ## SIG A, \
                       ElConstDistMatrix_ ## SIG x, \
    CREFLECT(F) beta,  ElDistMatrix_ ## SIG y ) \
  { EL_TRY( \
      Hemv( CReflect(uplo), \
            CReflect(alpha), *CReflect(A), *CReflect(x), \
            CReflect(beta), *CReflect(y) ) ) } \
  /* Her */ \
  ElError ElHer_ ## SIG \
  ( ElUpperOrLower uplo, Base<F> alpha, \
    ElConstMatrix_ ## SIG x, ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      Her( CReflect(uplo), alpha, \
           *CReflect(x), *CReflect(A) ) ) } \
  ElError ElHerDist_ ## SIG \
  ( ElUpperOrLower uplo, Base<F> alpha, \
    ElConstDistMatrix_ ## SIG x, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      Her( CReflect(uplo), alpha, *CReflect(x), *CReflect(A) ) ) } \
  /* Her2 */ \
  ElError ElHer2_ ## SIG \
  ( ElUpperOrLower uplo, CREFLECT(F) alpha, \
    ElConstMatrix_ ## SIG x, ElConstMatrix_ ## SIG y, \
    ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      Her2( CReflect(uplo), CReflect(alpha), \
            *CReflect(x), *CReflect(y), *CReflect(A) ) ) } \
  ElError ElHer2Dist_ ## SIG \
  ( ElUpperOrLower uplo, CREFLECT(F) alpha, \
    ElConstDistMatrix_ ## SIG x, ElConstDistMatrix_ ## SIG y, \
    ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      Her2( CReflect(uplo), CReflect(alpha), \
            *CReflect(x), *CReflect(y), *CReflect(A) ) ) } \
  /* Trr */ \
  ElError ElTrr_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(F) alpha, ElConstMatrix_ ## SIG x, ElConstMatrix_ ## SIG y, \
    ElMatrix_ ## SIG A, bool conjugate ) \
  { EL_TRY( \
      Trr( CReflect(uplo), \
           CReflect(alpha), *CReflect(x), *CReflect(y), \
           *CReflect(A), conjugate ) ) } \
  ElError ElTrrDist_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(F) alpha, ElConstDistMatrix_ ## SIG x, \
                       ElConstDistMatrix_ ## SIG y, \
                       ElDistMatrix_ ## SIG A, bool conjugate ) \
  { EL_TRY( \
      Trr( CReflect(uplo), \
           CReflect(alpha), *CReflect(x), *CReflect(y), *CReflect(A), \
           conjugate ) ) } \
  /* Trr2 */ \
  ElError ElTrr2_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(F) alpha, ElConstMatrix_ ## SIG X, ElConstMatrix_ ## SIG Y, \
    ElMatrix_ ## SIG A, bool conjugate ) \
  { EL_TRY( \
      Trr2( CReflect(uplo), \
            CReflect(alpha), *CReflect(X), *CReflect(Y), \
            *CReflect(A), conjugate ) ) } \
  ElError ElTrr2Dist_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(F) alpha, ElConstDistMatrix_ ## SIG X, \
                       ElConstDistMatrix_ ## SIG Y, \
                       ElDistMatrix_ ## SIG A, bool conjugate ) \
  { EL_TRY( \
      Trr2( CReflect(uplo), \
            CReflect(alpha), *CReflect(X), *CReflect(Y), *CReflect(A), \
            conjugate ) ) }

#include "El/macros/CInstantiate.h"

} // extern "C"
