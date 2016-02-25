/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El/core/FlamePart.hpp"

#include "El.h"
using namespace El;

extern "C" {

#define C_PROTO(SIG,SIGBASE,T) \
  /* Downward */ \
  ElError ElRepartitionDown_ ## SIG \
  ( ElMatrix_ ## SIG AT, ElMatrix_ ## SIG A0, \
                         ElMatrix_ ## SIG A1, \
    ElMatrix_ ## SIG AB, ElMatrix_ ## SIG A2, ElInt bsize ); \
  ElError ElRepartitionDownDist_ ## SIG \
  ( ElDistMatrix_ ## SIG AT, ElDistMatrix_ ## SIG A0, \
                             ElDistMatrix_ ## SIG A1, \
    ElDistMatrix_ ## SIG AB, ElDistMatrix_ ## SIG A2, ElInt bsize ); \
  ElError ElLockedRepartitionDown_ ## SIG \
  ( ElConstMatrix_ ## SIG AT, ElMatrix_ ## SIG A0, \
                              ElMatrix_ ## SIG A1, \
    ElConstMatrix_ ## SIG AB, ElMatrix_ ## SIG A2, ElInt bsize ); \
  ElError ElLockedRepartitionDownDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AT, ElDistMatrix_ ## SIG A0, \
                                  ElDistMatrix_ ## SIG A1, \
    ElConstDistMatrix_ ## SIG AB, ElDistMatrix_ ## SIG A2, ElInt bsize ); \
  /* Upward */ \
  ElError ElRepartitionUp_ ## SIG \
  ( ElMatrix_ ## SIG AT, ElMatrix_ ## SIG A0, \
                         ElMatrix_ ## SIG A1, \
    ElMatrix_ ## SIG AB, ElMatrix_ ## SIG A2, ElInt bsize ); \
  ElError ElRepartitionUpDist_ ## SIG \
  ( ElDistMatrix_ ## SIG AT, ElDistMatrix_ ## SIG A0, \
                             ElDistMatrix_ ## SIG A1, \
    ElDistMatrix_ ## SIG AB, ElDistMatrix_ ## SIG A2, ElInt bsize ); \
  ElError ElLockedRepartitionUp_ ## SIG \
  ( ElConstMatrix_ ## SIG AT, ElMatrix_ ## SIG A0, \
                              ElMatrix_ ## SIG A1, \
    ElConstMatrix_ ## SIG AB, ElMatrix_ ## SIG A2, ElInt bsize ); \
  ElError ElLockedRepartitionUpDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AT, ElDistMatrix_ ## SIG A0, \
                                  ElDistMatrix_ ## SIG A1, \
    ElConstDistMatrix_ ## SIG AB, ElDistMatrix_ ## SIG A2, ElInt bsize ); \
  /* Rightward */ \
  ElError ElRepartitionRight_ ## SIG \
  ( ElMatrix_ ## SIG AL, ElMatrix_ ## SIG AR, \
    ElMatrix_ ## SIG A0, ElMatrix_ ## SIG A1, ElMatrix_ ## SIG A2, \
    ElInt bsize ); \
  ElError ElRepartitionRightDist_ ## SIG \
  ( ElDistMatrix_ ## SIG AL, ElDistMatrix_ ## SIG AR, \
    ElDistMatrix_ ## SIG A0, ElDistMatrix_ ## SIG A1, ElDistMatrix_ ## SIG A2, \
    ElInt bsize ); \
  ElError ElLockedRepartitionRight_ ## SIG \
  ( ElConstMatrix_ ## SIG AL, ElConstMatrix_ ## SIG AR, \
    ElMatrix_ ## SIG A0, ElMatrix_ ## SIG A1, ElMatrix_ ## SIG A2, \
    ElInt bsize ); \
  ElError ElLockedRepartitionRightDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AL, ElConstDistMatrix_ ## SIG AR, \
    ElDistMatrix_ ## SIG A0, ElDistMatrix_ ## SIG A1, ElDistMatrix_ ## SIG A2, \
    ElInt bsize ); \
  /* Leftward */ \
  ElError ElRepartitionLeft_ ## SIG \
  ( ElMatrix_ ## SIG AL, ElMatrix_ ## SIG AR, \
    ElMatrix_ ## SIG A0, ElMatrix_ ## SIG A1, ElMatrix_ ## SIG A2, \
    ElInt bsize ); \
  ElError ElRepartitionLeftDist_ ## SIG \
  ( ElDistMatrix_ ## SIG AL, ElDistMatrix_ ## SIG AR, \
    ElDistMatrix_ ## SIG A0, ElDistMatrix_ ## SIG A1, ElDistMatrix_ ## SIG A2, \
    ElInt bsize ); \
  ElError ElLockedRepartitionLeft_ ## SIG \
  ( ElConstMatrix_ ## SIG AL, ElConstMatrix_ ## SIG AR, \
    ElMatrix_ ## SIG A0, ElMatrix_ ## SIG A1, ElMatrix_ ## SIG A2, \
    ElInt bsize ); \
  ElError ElLockedRepartitionLeftDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG AL, ElConstDistMatrix_ ## SIG AR, \
    ElDistMatrix_ ## SIG A0, ElDistMatrix_ ## SIG A1, ElDistMatrix_ ## SIG A2, \
    ElInt bsize ); \
  /* Down a diagonal */ \
  ElError ElRepartitionDownDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG ATL, ElMatrix_ ## SIG ATR, \
    ElMatrix_ ## SIG A00, ElMatrix_ ## SIG A01, ElMatrix_ ## SIG A02, \
    ElMatrix_ ## SIG A10, ElMatrix_ ## SIG A11, ElMatrix_ ## SIG A12, \
    ElMatrix_ ## SIG ABL, ElMatrix_ ## SIG ABR, \
    ElMatrix_ ## SIG A20, ElMatrix_ ## SIG A21, ElMatrix_ ## SIG A22, \
    ElInt bsize ); \
  ElError ElRepartitionDownDiagonalDist_ ## SIG \
  ( ElDistMatrix_ ## SIG ATL, ElDistMatrix_ ## SIG ATR, \
    ElDistMatrix_ ## SIG A00, ElDistMatrix_ ## SIG A01, ElDistMatrix_ ## SIG A02, \
    ElDistMatrix_ ## SIG A10, ElDistMatrix_ ## SIG A11, ElDistMatrix_ ## SIG A12, \
    ElDistMatrix_ ## SIG ABL, ElDistMatrix_ ## SIG ABR, \
    ElDistMatrix_ ## SIG A20, ElDistMatrix_ ## SIG A21, ElDistMatrix_ ## SIG A22, \
    ElInt bsize ); \
  ElError ElLockedRepartitionDownDiagonal_ ## SIG \
  ( ElConstMatrix_ ## SIG ATL, ElConstMatrix_ ## SIG ATR, \
    ElMatrix_ ## SIG A00, ElMatrix_ ## SIG A01, ElMatrix_ ## SIG A02, \
    ElMatrix_ ## SIG A10, ElMatrix_ ## SIG A11, ElMatrix_ ## SIG A12, \
    ElConstMatrix_ ## SIG ABL, ElConstMatrix_ ## SIG ABR, \
    ElMatrix_ ## SIG A20, ElMatrix_ ## SIG A21, ElMatrix_ ## SIG A22, \
    ElInt bsize ); \
  ElError ElLockedRepartitionDownDiagonalDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG ATL, ElConstDistMatrix_ ## SIG ATR, \
    ElDistMatrix_ ## SIG A00, ElDistMatrix_ ## SIG A01, ElDistMatrix_ ## SIG A02, \
    ElDistMatrix_ ## SIG A10, ElDistMatrix_ ## SIG A11, ElDistMatrix_ ## SIG A12, \
    ElConstDistMatrix_ ## SIG ABL, ElConstDistMatrix_ ## SIG ABR, \
    ElDistMatrix_ ## SIG A20, ElDistMatrix_ ## SIG A21, ElDistMatrix_ ## SIG A22, \
    ElInt bsize ); \
  /* Up a diagonal */ \
  ElError ElRepartitionUpDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG ATL, ElMatrix_ ## SIG ATR, \
    ElMatrix_ ## SIG A00, ElMatrix_ ## SIG A01, ElMatrix_ ## SIG A02, \
    ElMatrix_ ## SIG A10, ElMatrix_ ## SIG A11, ElMatrix_ ## SIG A12, \
    ElMatrix_ ## SIG ABL, ElMatrix_ ## SIG ABR, \
    ElMatrix_ ## SIG A20, ElMatrix_ ## SIG A21, ElMatrix_ ## SIG A22, \
    ElInt bsize ); \
  ElError ElRepartitionUpDiagonalDist_ ## SIG \
  ( ElDistMatrix_ ## SIG ATL, ElDistMatrix_ ## SIG ATR, \
    ElDistMatrix_ ## SIG A00, ElDistMatrix_ ## SIG A01, ElDistMatrix_ ## SIG A02, \
    ElDistMatrix_ ## SIG A10, ElDistMatrix_ ## SIG A11, ElDistMatrix_ ## SIG A12, \
    ElDistMatrix_ ## SIG ABL, ElDistMatrix_ ## SIG ABR, \
    ElDistMatrix_ ## SIG A20, ElDistMatrix_ ## SIG A21, ElDistMatrix_ ## SIG A22, \
    ElInt bsize ); \
  ElError ElLockedRepartitionUpDiagonal_ ## SIG \
  ( ElConstMatrix_ ## SIG ATL, ElConstMatrix_ ## SIG ATR, \
    ElMatrix_ ## SIG A00, ElMatrix_ ## SIG A01, ElMatrix_ ## SIG A02, \
    ElMatrix_ ## SIG A10, ElMatrix_ ## SIG A11, ElMatrix_ ## SIG A12, \
    ElConstMatrix_ ## SIG ABL, ElConstMatrix_ ## SIG ABR, \
    ElMatrix_ ## SIG A20, ElMatrix_ ## SIG A21, ElMatrix_ ## SIG A22, \
    ElInt bsize ); \
  ElError ElLockedRepartitionUpDiagonalDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG ATL, ElConstDistMatrix_ ## SIG ATR, \
    ElDistMatrix_ ## SIG A00, ElDistMatrix_ ## SIG A01, ElDistMatrix_ ## SIG A02, \
    ElDistMatrix_ ## SIG A10, ElDistMatrix_ ## SIG A11, ElDistMatrix_ ## SIG A12, \
    ElConstDistMatrix_ ## SIG ABL, ElConstDistMatrix_ ## SIG ABR, \
    ElDistMatrix_ ## SIG A20, ElDistMatrix_ ## SIG A21, ElDistMatrix_ ## SIG A22, \
    ElInt bsize );

#include "El/macros/CInstantiate.h"

} // extern "C"
