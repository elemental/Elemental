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

extern "C" {

#define C_PROTO_BASE(SIG,T) \
  /* Ones */ \
  ElError ElOnes_ ## SIG ( ElMatrix_ ## SIG A, ElInt m, ElInt n ) \
  { EL_TRY( Ones( *Reinterpret(A), m, n ) ) } \
  ElError ElOnesDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt m, ElInt n ) \
  { EL_TRY( Ones( *Reinterpret(A), m, n ) ) } \
  /* Uniform */ \
  ElError ElUniform_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt m, ElInt n, \
    CREFLECT(T) center, Base<T> radius ) \
  { EL_TRY( Uniform( *Reinterpret(A), m, n, Reinterpret(center), radius ) ) } \
  ElError ElUniformDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt m, ElInt n, \
    CREFLECT(T) center, Base<T> radius ) \
  { EL_TRY( Uniform( *Reinterpret(A), m, n, Reinterpret(center), radius ) ) }

#define C_PROTO_NOINT(SIG,T) \
  /* Cauchy */ \
  ElError ElCauchy_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt xSize, CREFLECT(T)* xBuf, \
                        ElInt ySize, CREFLECT(T)* yBuf ) \
  { try { \
      std::vector<T> x( Reinterpret(xBuf), Reinterpret(xBuf)+xSize ); \
      std::vector<T> y( Reinterpret(yBuf), Reinterpret(yBuf)+ySize ); \
      Cauchy( *Reinterpret(A), x, y ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElCauchyDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt xSize, CREFLECT(T)* xBuf, \
                            ElInt ySize, CREFLECT(T)* yBuf ) \
  { try { \
      std::vector<T> x( Reinterpret(xBuf), Reinterpret(xBuf)+xSize ); \
      std::vector<T> y( Reinterpret(yBuf), Reinterpret(yBuf)+ySize ); \
      Cauchy( *Reinterpret(A), x, y ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* Cauchy-like */ \
  ElError ElCauchyLike_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt rSize, CREFLECT(T)* rBuf, \
                        ElInt sSize, CREFLECT(T)* sBuf, \
                        ElInt xSize, CREFLECT(T)* xBuf, \
                        ElInt ySize, CREFLECT(T)* yBuf ) \
  { try { \
      std::vector<T> r( Reinterpret(rBuf), Reinterpret(rBuf)+rSize ); \
      std::vector<T> s( Reinterpret(sBuf), Reinterpret(sBuf)+sSize ); \
      std::vector<T> x( Reinterpret(xBuf), Reinterpret(xBuf)+xSize ); \
      std::vector<T> y( Reinterpret(yBuf), Reinterpret(yBuf)+ySize ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElCauchyLikeDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt rSize, CREFLECT(T)* rBuf, \
                            ElInt sSize, CREFLECT(T)* sBuf, \
                            ElInt xSize, CREFLECT(T)* xBuf, \
                            ElInt ySize, CREFLECT(T)* yBuf ) \
  { try { \
      std::vector<T> r( Reinterpret(rBuf), Reinterpret(rBuf)+rSize ); \
      std::vector<T> s( Reinterpret(sBuf), Reinterpret(sBuf)+sSize ); \
      std::vector<T> x( Reinterpret(xBuf), Reinterpret(xBuf)+xSize ); \
      std::vector<T> y( Reinterpret(yBuf), Reinterpret(yBuf)+ySize ); \
      CauchyLike( *Reinterpret(A), r, s, x, y ); \
    } EL_CATCH; return EL_SUCCESS; }

#define C_PROTO_INT(SIG,T) \
  C_PROTO_BASE(SIG,T)

#define C_PROTO_REAL(SIG,T) \
  C_PROTO_BASE(SIG,T) \
  C_PROTO_NOINT(SIG,T)

#define C_PROTO_COMPLEX(SIG,SIGBASE,T) \
  C_PROTO_BASE(SIG,T) \
  C_PROTO_NOINT(SIG,T) \
  /* Bull's head */ \
  ElError ElBullsHead_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( BullsHead( *Reinterpret(A), n ) ) } \
  ElError ElBullsHeadDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( BullsHead( *Reinterpret(A), n ) ) } 

#include "El/macros/CInstantiate.h"

} // extern "C"
