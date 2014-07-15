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

#define C_PROTO(SIG,T) \
  /* B = A */ \
  ElError ElCopyDistMatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElDistMatrix_ ## SIG BHandle ) \
  { EL_TRY( Copy( *Reinterpret(AHandle), *Reinterpret(BHandle) ) ) } \
  /* B = A^T */ \
  ElError ElTransposeDistMatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElDistMatrix_ ## SIG BHandle ) \
  { EL_TRY( Transpose(*Reinterpret(AHandle),*Reinterpret(BHandle),false) ) } \
  /* B = A^H */ \
  ElError ElAdjointDistMatrix_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AHandle, ElDistMatrix_ ## SIG BHandle ) \
  { EL_TRY( Adjoint( *Reinterpret(AHandle), *Reinterpret(BHandle) ) ) }

#include "El/macros/CInstantiate.h"

} // extern "C"
