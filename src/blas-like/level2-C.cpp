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
            Reinterpret(beta), DM_CAST(T,y) ) ) } \
  /* TODO: Trmv */ \
  /* TODO: Trr2 */ \
  /* TODO: Trr */ 

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

/* TODO: QuasiTrsv */
/* TODO: Trsv */

#define C_PROTO_INT(SIG,T) C_PROTO_BASE(SIG,T)

#define C_PROTO(SIG,T) \
  C_PROTO_BASE(SIG,T) \
  C_PROTO_NOINT(SIG,T)

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
           DM_CAST_CONST(T,x), DM_CAST(T,A) ) ) }

#include "El/macros/CInstantiate.h"

} // extern "C"
