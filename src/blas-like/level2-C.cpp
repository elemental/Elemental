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

#define M_CAST(T,A) dynamic_cast<Matrix<T>&>(*Reinterpret(A))
#define M_CAST_CONST(T,A) dynamic_cast<const Matrix<T>&>(*Reinterpret(A))

#define DM_CAST(T,A) dynamic_cast<DistMatrix<T>&>(*Reinterpret(A))
#define DM_CAST_CONST(T,A) dynamic_cast<const DistMatrix<T>&>(*Reinterpret(A))

extern "C" {

#define C_PROTO_BASE(SIG,T) \
  /* Gemv */ \
  ElError ElGemv_ ## SIG \
  ( ElOrientation orientation, \
    CREFLECT(T) alpha, ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG x, \
    CREFLECT(T) beta, ElMatrix_ ## SIG y ) \
  { EL_TRY( \
      Gemv( Reinterpret(orientation), \
            Reinterpret(alpha), M_CAST_CONST(T,A), M_CAST_CONST(T,x), \
            Reinterpret(beta), M_CAST(T,y) ) ) } \
  ElError ElGemvDist_ ## SIG \
  ( ElOrientation orientation, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG A, \
                       ElConstDistMatrix_ ## SIG x, \
    CREFLECT(T) beta,  ElDistMatrix_ ## SIG y ) \
  { EL_TRY( \
      Gemv( Reinterpret(orientation), \
            Reinterpret(alpha), DM_CAST_CONST(T,A), DM_CAST_CONST(T,x), \
            Reinterpret(beta), DM_CAST(T,y) ) ) }

#define C_PROTO_NOINT(SIG,T) \
  /* Ger */ \
  ElError ElGer_ ## SIG \
  ( CREFLECT(T) alpha, ElConstMatrix_ ## SIG x, ElConstMatrix_ ## SIG y, \
    ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      Ger( Reinterpret(alpha), M_CAST_CONST(T,x), M_CAST_CONST(T,y), \
           M_CAST(T,A) ) ) } \
  ElError ElGerDist_ ## SIG \
  ( CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG x, \
                       ElConstDistMatrix_ ## SIG y, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      Ger( Reinterpret(alpha), DM_CAST_CONST(T,x), DM_CAST_CONST(T,y), \
           DM_CAST(T,A) ) ) } \
  /* Geru */ \
  ElError ElGeru_ ## SIG \
  ( CREFLECT(T) alpha, ElConstMatrix_ ## SIG x, ElConstMatrix_ ## SIG y, \
    ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      Geru( Reinterpret(alpha), M_CAST_CONST(T,x), M_CAST_CONST(T,y), \
           M_CAST(T,A) ) ) } \
  ElError ElGeruDist_ ## SIG \
  ( CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG x, \
                       ElConstDistMatrix_ ## SIG y, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      Geru( Reinterpret(alpha), DM_CAST_CONST(T,x), DM_CAST_CONST(T,y), \
            DM_CAST(T,A) ) ) } \
  /* Symv */ \
  ElError ElSymv_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG x, \
    CREFLECT(T) beta, ElMatrix_ ## SIG y ) \
  { EL_TRY( \
      Symv( Reinterpret(uplo), \
            Reinterpret(alpha), M_CAST_CONST(T,A), M_CAST_CONST(T,x), \
            Reinterpret(beta), M_CAST(T,y) ) ) } \
  ElError ElSymvDist_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG A, \
                       ElConstDistMatrix_ ## SIG x, \
    CREFLECT(T) beta,  ElDistMatrix_ ## SIG y ) \
  { EL_TRY( \
      Symv( Reinterpret(uplo), \
            Reinterpret(alpha), DM_CAST_CONST(T,A), DM_CAST_CONST(T,x), \
            Reinterpret(beta), DM_CAST(T,y) ) ) } \
  /* Syr */ \
  ElError ElSyr_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstMatrix_ ## SIG x, ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      Syr( Reinterpret(uplo), \
            Reinterpret(alpha), M_CAST_CONST(T,x), M_CAST(T,A) ) ) } \
  ElError ElSyrDist_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG x, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      Syr( Reinterpret(uplo), \
            Reinterpret(alpha), DM_CAST_CONST(T,x), DM_CAST(T,A) ) ) } \
  /* Syr2 */ \
  ElError ElSyr2_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstMatrix_ ## SIG x, ElConstMatrix_ ## SIG y, \
    ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      Syr2( Reinterpret(uplo), \
            Reinterpret(alpha), M_CAST_CONST(T,x), M_CAST_CONST(T,y), \
            M_CAST(T,A) ) ) } \
  ElError ElSyr2Dist_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG x, \
                       ElConstDistMatrix_ ## SIG y, \
                       ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      Syr2( Reinterpret(uplo), \
            Reinterpret(alpha), DM_CAST_CONST(T,x), DM_CAST_CONST(T,y), \
            DM_CAST(T,A) ) ) } \
  /* QuasiTrsv */ \
  ElError ElQuasiTrsv_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, ElConstMatrix_ ## SIG A, \
    ElMatrix_ ## SIG x ) \
  { EL_TRY( \
      QuasiTrsv( \
        Reinterpret(uplo), Reinterpret(orientation), \
        M_CAST_CONST(T,A), M_CAST(T,x) ) ) } \
  ElError ElQuasiTrsvDist_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG x ) \
  { EL_TRY( \
      QuasiTrsv( \
        Reinterpret(uplo), Reinterpret(orientation), \
        M_CAST_CONST(T,A), M_CAST(T,x) ) ) } \
  /* Trmv */ \
  ElError ElTrmv_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, ElUnitOrNonUnit diag, \
    ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG x ) \
  { EL_TRY( \
      Trmv( \
        Reinterpret(uplo), Reinterpret(orientation), Reinterpret(diag), \
        M_CAST_CONST(T,A), M_CAST(T,x) ) ) } \
  /* NOTE: Distributed Trmv not implemented, use Trmm */ \
  /* Trsv */ \
  ElError ElTrsv_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, ElUnitOrNonUnit diag, \
    ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG x ) \
  { EL_TRY( \
      Trsv( \
        Reinterpret(uplo), Reinterpret(orientation), Reinterpret(diag), \
        M_CAST_CONST(T,A), M_CAST(T,x) ) ) } \
  ElError ElTrsvDist_ ## SIG \
  ( ElUpperOrLower uplo, ElOrientation orientation, ElUnitOrNonUnit diag, \
    ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG x ) \
  { EL_TRY( \
      Trsv( \
        Reinterpret(uplo), Reinterpret(orientation), Reinterpret(diag), \
        M_CAST_CONST(T,A), M_CAST(T,x) ) ) }

#define C_PROTO_INT(SIG,T) C_PROTO_BASE(SIG,T)

#define C_PROTO(SIG,T) \
  C_PROTO_BASE(SIG,T) \
  C_PROTO_NOINT(SIG,T) \
  /* Trr */ \
  ElError ElTrr_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstMatrix_ ## SIG x, ElConstMatrix_ ## SIG y, \
    ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      Trr( Reinterpret(uplo), \
           Reinterpret(alpha), M_CAST_CONST(T,x), M_CAST_CONST(T,y), \
           M_CAST(T,A) ) ) } \
  ElError ElTrrDist_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG x, \
                       ElConstDistMatrix_ ## SIG y, \
                       ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      Trr( Reinterpret(uplo), \
           Reinterpret(alpha), DM_CAST_CONST(T,x), DM_CAST_CONST(T,y), \
           DM_CAST(T,A) ) ) } \
  /* Trr2 */ \
  ElError ElTrr2_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstMatrix_ ## SIG X, ElConstMatrix_ ## SIG Y, \
    ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      Trr2( Reinterpret(uplo), \
            Reinterpret(alpha), M_CAST_CONST(T,X), M_CAST_CONST(T,Y), \
            M_CAST(T,A) ) ) } \
  ElError ElTrr2Dist_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG X, \
                       ElConstDistMatrix_ ## SIG Y, \
                       ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      Trr2( Reinterpret(uplo), \
            Reinterpret(alpha), DM_CAST_CONST(T,X), DM_CAST_CONST(T,Y), \
            DM_CAST(T,A) ) ) } 

#define C_PROTO_COMPLEX(SIG,SIGBASE,T) \
  C_PROTO_BASE(SIG,T) \
  C_PROTO_NOINT(SIG,T) \
  /* Hemv */ \
  ElError ElHemv_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG x, \
    CREFLECT(T) beta, ElMatrix_ ## SIG y ) \
  { EL_TRY( \
      Hemv( Reinterpret(uplo), \
            Reinterpret(alpha), M_CAST_CONST(T,A), M_CAST_CONST(T,x), \
            Reinterpret(beta), M_CAST(T,y) ) ) } \
  ElError ElHemvDist_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG A, \
                       ElConstDistMatrix_ ## SIG x, \
    CREFLECT(T) beta,  ElDistMatrix_ ## SIG y ) \
  { EL_TRY( \
      Hemv( Reinterpret(uplo), \
            Reinterpret(alpha), DM_CAST_CONST(T,A), DM_CAST_CONST(T,x), \
            Reinterpret(beta), DM_CAST(T,y) ) ) } \
  /* Her */ \
  ElError ElHer_ ## SIG \
  ( ElUpperOrLower uplo, CREFLECT(T) alpha, \
    ElConstMatrix_ ## SIG x, ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      Her( Reinterpret(uplo), Reinterpret(alpha), \
           M_CAST_CONST(T,x), M_CAST(T,A) ) ) } \
  ElError ElHerDist_ ## SIG \
  ( ElUpperOrLower uplo, CREFLECT(T) alpha, \
    ElConstDistMatrix_ ## SIG x, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      Her( Reinterpret(uplo), Reinterpret(alpha), \
           DM_CAST_CONST(T,x), DM_CAST(T,A) ) ) } \
  /* Her2 */ \
  ElError ElHer2_ ## SIG \
  ( ElUpperOrLower uplo, CREFLECT(T) alpha, \
    ElConstMatrix_ ## SIG x, ElConstMatrix_ ## SIG y, \
    ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      Her2( Reinterpret(uplo), Reinterpret(alpha), \
            M_CAST_CONST(T,x), M_CAST_CONST(T,y), M_CAST(T,A) ) ) } \
  ElError ElHer2Dist_ ## SIG \
  ( ElUpperOrLower uplo, CREFLECT(T) alpha, \
    ElConstDistMatrix_ ## SIG x, ElConstDistMatrix_ ## SIG y, \
    ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      Her2( Reinterpret(uplo), Reinterpret(alpha), \
            DM_CAST_CONST(T,x), DM_CAST_CONST(T,y), DM_CAST(T,A) ) ) } \
  /* Trr */ \
  ElError ElTrr_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstMatrix_ ## SIG x, ElConstMatrix_ ## SIG y, \
    ElMatrix_ ## SIG A, bool conjugate ) \
  { EL_TRY( \
      Trr( Reinterpret(uplo), \
           Reinterpret(alpha), M_CAST_CONST(T,x), M_CAST_CONST(T,y), \
           M_CAST(T,A), conjugate ) ) } \
  ElError ElTrrDist_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG x, \
                       ElConstDistMatrix_ ## SIG y, \
                       ElDistMatrix_ ## SIG A, bool conjugate ) \
  { EL_TRY( \
      Trr( Reinterpret(uplo), \
           Reinterpret(alpha), DM_CAST_CONST(T,x), DM_CAST_CONST(T,y), \
           DM_CAST(T,A), conjugate ) ) } \
  /* Trr2 */ \
  ElError ElTrr2_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstMatrix_ ## SIG X, ElConstMatrix_ ## SIG Y, \
    ElMatrix_ ## SIG A, bool conjugate ) \
  { EL_TRY( \
      Trr2( Reinterpret(uplo), \
            Reinterpret(alpha), M_CAST_CONST(T,X), M_CAST_CONST(T,Y), \
            M_CAST(T,A), conjugate ) ) } \
  ElError ElTrr2Dist_ ## SIG \
  ( ElUpperOrLower uplo, \
    CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG X, \
                       ElConstDistMatrix_ ## SIG Y, \
                       ElDistMatrix_ ## SIG A, bool conjugate ) \
  { EL_TRY( \
      Trr2( Reinterpret(uplo), \
            Reinterpret(alpha), DM_CAST_CONST(T,X), DM_CAST_CONST(T,Y), \
            DM_CAST(T,A), conjugate ) ) }

#include "El/macros/CInstantiate.h"

} // extern "C"
