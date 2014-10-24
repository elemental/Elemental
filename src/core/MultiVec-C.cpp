/*
   Copyright (c) 2009-2014, Jack Poulson
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
  /* Constructors and destructors */ \
  ElError ElMultiVecCreate_ ## SIG ( ElMultiVec_ ## SIG * A ) \
  { EL_TRY( *A = CReflect( new MultiVec<T>() ) ) } \
  ElError ElMultiVecDestroy_ ## SIG ( ElConstMultiVec_ ## SIG A ) \
  { EL_TRY( delete CReflect(A) ) } \
  /* Assignment and reconfiguration */ \
  ElError ElMultiVecEmpty_ ## SIG ( ElMultiVec_ ## SIG A ) \
  { EL_TRY( CReflect(A)->Empty() ) } \
  ElError ElMultiVecResize_ ## SIG \
  ( ElMultiVec_ ## SIG A, ElInt height, ElInt width ) \
  { EL_TRY( CReflect(A)->Resize(height,width) ) } \
  /* Queries */ \
  ElError ElMultiVecHeight_ ## SIG \
  ( ElConstMultiVec_ ## SIG A, ElInt* height ) \
  { EL_TRY( *height = CReflect(A)->Height() ) } \
  ElError ElMultiVecWidth_ ## SIG \
  ( ElConstMultiVec_ ## SIG A, ElInt* width ) \
  { EL_TRY( *width = CReflect(A)->Width() ) } \
  /* Entrywise manipulation */ \
  ElError ElMultiVecGet_ ## SIG \
  ( ElMultiVec_ ## SIG A, ElInt i, ElInt j, CREFLECT(T)* value ) \
  { EL_TRY( *value = CReflect(CReflect(A)->Get(i,j)) ) } \
  ElError ElMultiVecSet_ ## SIG \
  ( ElMultiVec_ ## SIG A, ElInt i, ElInt j, CREFLECT(T) value ) \
  { EL_TRY( CReflect(A)->Set(i,j,CReflect(value)) ) } \
  ElError ElMultiVecUpdate_ ## SIG \
  ( ElMultiVec_ ## SIG A, ElInt i, ElInt j, CREFLECT(T) value ) \
  { EL_TRY( CReflect(A)->Update(i,j,CReflect(value)) ) }

#include "El/macros/CInstantiate.h"

} // extern "C"
