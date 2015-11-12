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

#define C_PROTO(SIG,SIGBASE,F) \
  ElError ElLatticeGramSchmidt_ ## SIG \
  ( ElConstMatrix_ ## SIG B, ElMatrix_ ## SIG G, ElMatrix_ ## SIG M ) \
  { EL_TRY( LatticeGramSchmidt( *CReflect(B), *CReflect(G), *CReflect(M) ) ) } \
  ElError ElLLL_ ## SIG \
  ( ElMatrix_ ## SIG B, \
    ElMatrix_ ## SIG QR, \
    Base<F> delta, \
    Base<F> innerTol, \
    bool presort, \
    bool smallestFirst, \
    bool progress, \
    ElInt* numBacktrack ) \
  { EL_TRY( *numBacktrack = \
      LLL( *CReflect(B), *CReflect(QR), delta, innerTol, \
           presort, smallestFirst, progress ) ) } \
  ElError ElLLLDelta_ ## SIG \
  ( ElConstMatrix_ ## SIG QR, Base<F>* delta ) \
  { EL_TRY( *delta = LLLDelta( *CReflect(QR) ) ) }

#define EL_NO_INT_PROTO
#include "El/macros/CInstantiate.h"

} // extern "C"
