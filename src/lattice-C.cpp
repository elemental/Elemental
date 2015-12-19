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

ElError ElLLLCtrlDefault_s( ElLLLCtrl_s* ctrl )
{
    ctrl->delta = 0.75f;
    ctrl->weak = false;
    ctrl->presort = true;
    ctrl->smallestFirst = true;
    ctrl->reorthogTol = 0;
    ctrl->zeroTol = limits::Epsilon<float>();
    ctrl->progress = false;
    ctrl->time = false;
    return EL_SUCCESS;
}

ElError ElLLLCtrlDefault_d( ElLLLCtrl_d* ctrl )
{
    ctrl->delta = 0.75;
    ctrl->weak = false;
    ctrl->presort = true;
    ctrl->smallestFirst = true;
    ctrl->reorthogTol = 0;
    ctrl->zeroTol = limits::Epsilon<double>();
    ctrl->progress = false;
    ctrl->time = false;
    return EL_SUCCESS;
}

#define C_PROTO(SIG,SIGBASE,F) \
  ElError ElLatticeGramSchmidt_ ## SIG \
  ( ElConstMatrix_ ## SIG B, ElMatrix_ ## SIG G, ElMatrix_ ## SIG M ) \
  { EL_TRY( LatticeGramSchmidt( *CReflect(B), *CReflect(G), *CReflect(M) ) ) } \
  ElError ElLLL_ ## SIG \
  ( ElMatrix_ ## SIG B, \
    ElMatrix_ ## SIG QR, \
    ElLLLCtrl_ ## SIGBASE ctrl, \
    ElLLLInfo* infoC ) \
  { EL_TRY( \
      auto info = LLL( *CReflect(B), *CReflect(QR), CReflect(ctrl) );\
      *infoC = CReflect(info); \
    ) } \
  ElError ElLLLFull_ ## SIG \
  ( ElMatrix_ ## SIG B, \
    ElMatrix_ ## SIG U, \
    ElMatrix_ ## SIG UInv, \
    ElMatrix_ ## SIG QR, \
    ElLLLCtrl_ ## SIGBASE ctrl, \
    ElLLLInfo* infoC ) \
  { EL_TRY( \
      auto info = \
        LLL( *CReflect(B), *CReflect(U), *CReflect(UInv), *CReflect(QR), \
             CReflect(ctrl) ); \
      *infoC = CReflect(info); \
    ) } \
  ElError ElLLLDelta_ ## SIG \
  ( ElConstMatrix_ ## SIG QR, ElLLLCtrl_ ## SIGBASE ctrl, Base<F>* delta ) \
  { EL_TRY( *delta = LLLDelta( *CReflect(QR), CReflect(ctrl) ) ) }

#define EL_NO_INT_PROTO
#include "El/macros/CInstantiate.h"

} // extern "C"
