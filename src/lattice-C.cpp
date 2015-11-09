/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El.h"
using namespace El;

extern "C" {

#define C_PROTO_REAL(SIG,Real) \
  ElError ElLLL_ ## SIG \
  ( ElMatrix_ ## SIG B, Real delta, Real eta, Real theta, Real innerTol ) \
  { EL_TRY( LLL( *CReflect(B), delta, eta, theta, innerTol ) ) }

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/CInstantiate.h"

} // extern "C"
