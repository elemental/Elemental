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

#define C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* Basic complex entry manipulation
     ================================ */ \
  ElError ElArg_ ## SIG ( CREFLECT(F) alpha, Base<F>* alphaArg ) \
  { EL_TRY( *alphaArg = Arg(CReflect(alpha)) ) } \
  /* Size measurement
     ================ */ \
  ElError ElAbs_ ## SIG ( CREFLECT(F) alpha, Base<F>* alphaAbs ) \
  { EL_TRY( *alphaAbs = Abs(CReflect(alpha)) ) } \
  /* Exponentiation 
     ============== */ \
  ElError ElExp_ ## SIG ( CREFLECT(F) alpha, CREFLECT(F)* alphaExp ) \
  { EL_TRY( *alphaExp = CReflect(Exp(CReflect(alpha))) ) } \
  ElError ElPow_ ## SIG \
  ( CREFLECT(F) alpha, CREFLECT(F) beta, CREFLECT(F)* result ) \
  { EL_TRY( *result = CReflect(Pow(CReflect(alpha),CReflect(beta))) ) } \
  /* Inverse exponentiation
     ---------------------- */ \
  ElError ElLog_ ## SIG ( CREFLECT(F) alpha, CREFLECT(F)* alphaLog ) \
  { EL_TRY( *alphaLog = CReflect(Log(CReflect(alpha))) ) } \
  ElError ElSqrt_ ## SIG ( CREFLECT(F) alpha, CREFLECT(F)* alphaSqrt ) \
  { EL_TRY( *alphaSqrt = CReflect(Sqrt(CReflect(alpha))) ) } \
  /* Trigonometric
     ============= */ \
  ElError ElCos_ ## SIG ( CREFLECT(F) alpha, CREFLECT(F)* alphaCos ) \
  { EL_TRY( *alphaCos = CReflect(Cos(CReflect(alpha))) ) } \
  ElError ElSin_ ## SIG ( CREFLECT(F) alpha, CREFLECT(F)* alphaSin ) \
  { EL_TRY( *alphaSin = CReflect(Sin(CReflect(alpha))) ) } \
  ElError ElTan_ ## SIG ( CREFLECT(F) alpha, CREFLECT(F)* alphaTan ) \
  { EL_TRY( *alphaTan = CReflect(Tan(CReflect(alpha))) ) } \
  /* Inverse trigonometric
     --------------------- */ \
  ElError ElAcos_ ## SIG ( CREFLECT(F) alpha, CREFLECT(F)* alphaAcos ) \
  { EL_TRY( *alphaAcos = CReflect(Acos(CReflect(alpha))) ) } \
  ElError ElAsin_ ## SIG ( CREFLECT(F) alpha, CREFLECT(F)* alphaAsin ) \
  { EL_TRY( *alphaAsin = CReflect(Asin(CReflect(alpha))) ) } \
  ElError ElAtan_ ## SIG ( CREFLECT(F) alpha, CREFLECT(F)* alphaAtan ) \
  { EL_TRY( *alphaAtan = CReflect(Atan(CReflect(alpha))) ) } \
  /* Hyperbolic 
     ========== */ \
  ElError ElCosh_ ## SIG ( CREFLECT(F) alpha, CREFLECT(F)* alphaCosh ) \
  { EL_TRY( *alphaCosh = CReflect(Cosh(CReflect(alpha))) ) } \
  ElError ElSinh_ ## SIG ( CREFLECT(F) alpha, CREFLECT(F)* alphaSinh ) \
  { EL_TRY( *alphaSinh = CReflect(Sinh(CReflect(alpha))) ) } \
  ElError ElTanh_ ## SIG ( CREFLECT(F) alpha, CREFLECT(F)* alphaTanh ) \
  { EL_TRY( *alphaTanh = CReflect(Tanh(CReflect(alpha))) ) } \
  /* Inverse hyperbolic
     ------------------ */ \
  ElError ElAcosh_ ## SIG ( CREFLECT(F) alpha, CREFLECT(F)* alphaAcosh ) \
  { EL_TRY( *alphaAcosh = CReflect(Acosh(CReflect(alpha))) ) } \
  ElError ElAsinh_ ## SIG ( CREFLECT(F) alpha, CREFLECT(F)* alphaAsinh ) \
  { EL_TRY( *alphaAsinh = CReflect(Asinh(CReflect(alpha))) ) } \
  ElError ElAtanh_ ## SIG ( CREFLECT(F) alpha, CREFLECT(F)* alphaAtanh ) \
  { EL_TRY( *alphaAtanh = CReflect(Atanh(CReflect(alpha))) ) }

#define C_PROTO_REAL_ONLY(SIG,Real) \
  ElError ElAtan2_ ## SIG ( Real y, Real x, Real* result ) \
  { EL_TRY( *result = Atan2( y, x ) ) } \
  ElError ElSgn_ ## SIG ( Real alpha, bool symmetric, Real* result ) \
  { EL_TRY( *result = Sgn(alpha,symmetric) ) }

#define C_PROTO_COMPLEX_ONLY(SIG,SIGBASE,F) \
  /* Basic complex entry manipulation
     ================================ */ \
  ElError ElComplexFromPolar_ ## SIG \
  ( Base<F> r, Base<F> theta, CREFLECT(F)* alpha ) \
  { EL_TRY( *alpha = CReflect(ComplexFromPolar(r,theta)) ) } \
  /* Size measurement
     ================ */ \
  ElError ElSafeAbs_ ## SIG ( CREFLECT(F) alpha, Base<F>* alphaAbs ) \
  { EL_TRY( *alphaAbs = SafeAbs(CReflect(alpha)) ) } \
  ElError ElFastAbs_ ## SIG ( CREFLECT(F) alpha, Base<F>* alphaAbs ) \
  { EL_TRY( *alphaAbs = FastAbs(CReflect(alpha)) ) }

#define C_PROTO_INT(SIG,T) \
  ElError ElAbs_i( T alpha, T* alphaAbs ) \
  { EL_TRY( *alphaAbs = Abs(alpha) ) } \
  ElError ElSgn_i( T alpha, bool symmetric, T* result ) \
  { EL_TRY( *result = Sgn(alpha,symmetric) ) }

#define C_PROTO_REAL(SIG,F) \
  C_PROTO_FIELD(SIG,SIG,F) \
  C_PROTO_REAL_ONLY(SIG,F)

#define C_PROTO_COMPLEX(SIG,SIGBASE,F) \
  C_PROTO_FIELD(SIG,SIGBASE,F) \
  C_PROTO_COMPLEX_ONLY(SIG,SIGBASE,F)

#include "El/macros/CInstantiate.h"

} // extern "C"
