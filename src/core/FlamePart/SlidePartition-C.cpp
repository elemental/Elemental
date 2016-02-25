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
  ElError ElSlidePartitionDown_ ## SIG \
  ( ElMatrix_ ## SIG AT, ElMatrix_ ## SIG A0, \
                         ElMatrix_ ## SIG A1, \
    ElMatrix_ ## SIG AB, ElMatrix_ ## SIG A2 ); \
  ElError ElSlidePartitionDownDist_ ## SIG \
  ( ElDistMatrix_ ## SIG AT, ElDistMatrix_ ## SIG A0, \
                             ElDistMatrix_ ## SIG A1, \
    ElDistMatrix_ ## SIG AB, ElDistMatrix_ ## SIG A2 ); \
  ElError ElSlideLockedPartitionDown_ ## SIG \
  ( ElMatrix_ ## SIG AT, ElConstMatrix_ ## SIG A0, \
                         ElConstMatrix_ ## SIG A1, \
    ElMatrix_ ## SIG AB, ElConstMatrix_ ## SIG A2 ); \
  ElError ElSlideLockedPartitionDownDist_ ## SIG \
  ( ElDistMatrix_ ## SIG AT, ElConstDistMatrix_ ## SIG A0, \
                             ElConstDistMatrix_ ## SIG A1, \
    ElDistMatrix_ ## SIG AB, ElConstDistMatrix_ ## SIG A2 ); \
  /* Upward */ \
  ElError ElSlidePartitionUp_ ## SIG \
  ( ElMatrix_ ## SIG AT, ElMatrix_ ## SIG A0, \
                         ElMatrix_ ## SIG A1, \
    ElMatrix_ ## SIG AB, ElMatrix_ ## SIG A2 ); \
  ElError ElSlidePartitionUpDist_ ## SIG \
  ( ElDistMatrix_ ## SIG AT, ElDistMatrix_ ## SIG A0, \
                             ElDistMatrix_ ## SIG A1, \
    ElDistMatrix_ ## SIG AB, ElDistMatrix_ ## SIG A2 ); \
  ElError ElSlideLockedPartitionUp_ ## SIG \
  ( ElMatrix_ ## SIG AT, ElConstMatrix_ ## SIG A0, \
                         ElConstMatrix_ ## SIG A1, \
    ElMatrix_ ## SIG AB, ElConstMatrix_ ## SIG A2 ); \
  ElError ElSlideLockedPartitionUpDist_ ## SIG \
  ( ElDistMatrix_ ## SIG AT, ElConstDistMatrix_ ## SIG A0, \
                             ElConstDistMatrix_ ## SIG A1, \
    ElDistMatrix_ ## SIG AB, ElConstDistMatrix_ ## SIG A2 ); \
  /* Rightward */ \
  ElError ElSlidePartitionRight_ ## SIG \
  ( ElMatrix_ ## SIG AL, ElMatrix_ ## SIG AR, \
    ElMatrix_ ## SIG A0, ElMatrix_ ## SIG A1, ElMatrix_ ## SIG A2 ); \
  ElError ElSlidePartitionRightDist_ ## SIG \
  ( ElDistMatrix_ ## SIG AL, ElDistMatrix_ ## SIG AR, \
    ElDistMatrix_ ## SIG A0, ElDistMatrix_ ## SIG A1, ElDistMatrix_ ## SIG A2 ); \
  ElError ElSlideLockedPartitionRight_ ## SIG \
  ( ElMatrix_ ## SIG AL, ElMatrix_ ## SIG AR, \
    ElConstMatrix_ ## SIG A0, ElConstMatrix_ ## SIG A1, ElConstMatrix_ ## SIG A2 ); \
  ElError ElSlideLockedPartitionRightDist_ ## SIG \
  ( ElDistMatrix_ ## SIG AL, ElDistMatrix_ ## SIG AR, \
    ElConstDistMatrix_ ## SIG A0, ElConstDistMatrix_ ## SIG A1, ElConstDistMatrix_ ## SIG A2 ); \
  /* Leftward */ \
  ElError ElSlidePartitionLeft_ ## SIG \
  ( ElMatrix_ ## SIG AL, ElMatrix_ ## SIG AR, \
    ElMatrix_ ## SIG A0, ElMatrix_ ## SIG A1, ElMatrix_ ## SIG A2 ); \
  ElError ElSlidePartitionLeftDist_ ## SIG \
  ( ElDistMatrix_ ## SIG AL, ElDistMatrix_ ## SIG AR, \
    ElDistMatrix_ ## SIG A0, ElDistMatrix_ ## SIG A1, ElDistMatrix_ ## SIG A2 ); \
  ElError ElSlideLockedPartitionLeft_ ## SIG \
  ( ElMatrix_ ## SIG AL, ElMatrix_ ## SIG AR, \
    ElConstMatrix_ ## SIG A0, ElConstMatrix_ ## SIG A1, ElConstMatrix_ ## SIG A2 ); \
  ElError ElSlideLockedPartitionLeftDist_ ## SIG \
  ( ElDistMatrix_ ## SIG AL, ElDistMatrix_ ## SIG AR, \
    ElConstDistMatrix_ ## SIG A0, ElConstDistMatrix_ ## SIG A1, ElConstDistMatrix_ ## SIG A2 ); \
  /* Down a diagonal */ \
  ElError ElSlidePartitionDownDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG ATL, ElMatrix_ ## SIG ATR, \
    ElMatrix_ ## SIG A00, ElMatrix_ ## SIG A01, ElMatrix_ ## SIG A02, \
    ElMatrix_ ## SIG A10, ElMatrix_ ## SIG A11, ElMatrix_ ## SIG A12, \
    ElMatrix_ ## SIG ABL, ElMatrix_ ## SIG ABR, \
    ElMatrix_ ## SIG A20, ElMatrix_ ## SIG A21, ElMatrix_ ## SIG A22 ); \
  ElError ElSlidePartitionDownDiagonalDist_ ## SIG \
  ( ElDistMatrix_ ## SIG ATL, ElDistMatrix_ ## SIG ATR, \
    ElDistMatrix_ ## SIG A00, ElDistMatrix_ ## SIG A01, ElDistMatrix_ ## SIG A02, \
    ElDistMatrix_ ## SIG A10, ElDistMatrix_ ## SIG A11, ElDistMatrix_ ## SIG A12, \
    ElDistMatrix_ ## SIG ABL, ElDistMatrix_ ## SIG ABR, \
    ElDistMatrix_ ## SIG A20, ElDistMatrix_ ## SIG A21, ElDistMatrix_ ## SIG A22 ); \
  ElError ElSlideLockedPartitionDownDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG ATL, ElMatrix_ ## SIG ATR, \
    ElConstMatrix_ ## SIG A00, ElConstMatrix_ ## SIG A01, ElConstMatrix_ ## SIG A02, \
    ElConstMatrix_ ## SIG A10, ElConstMatrix_ ## SIG A11, ElConstMatrix_ ## SIG A12, \
    ElMatrix_ ## SIG ABL, ElMatrix_ ## SIG ABR, \
    ElConstMatrix_ ## SIG A20, ElConstMatrix_ ## SIG A21, ElConstMatrix_ ## SIG A22 ); \
  ElError ElSlideLockedPartitionDownDiagonalDist_ ## SIG \
  ( ElDistMatrix_ ## SIG ATL, ElDistMatrix_ ## SIG ATR, \
    ElConstDistMatrix_ ## SIG A00, ElConstDistMatrix_ ## SIG A01, ElConstDistMatrix_ ## SIG A02, \
    ElConstDistMatrix_ ## SIG A10, ElConstDistMatrix_ ## SIG A11, ElConstDistMatrix_ ## SIG A12, \
    ElDistMatrix_ ## SIG ABL, ElDistMatrix_ ## SIG ABR, \
    ElConstDistMatrix_ ## SIG A20, ElConstDistMatrix_ ## SIG A21, ElConstDistMatrix_ ## SIG A22 ); \
  /* Up a diagonal */ \
  ElError ElSlidePartitionUpDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG ATL, ElMatrix_ ## SIG ATR, \
    ElMatrix_ ## SIG A00, ElMatrix_ ## SIG A01, ElMatrix_ ## SIG A02, \
    ElMatrix_ ## SIG A10, ElMatrix_ ## SIG A11, ElMatrix_ ## SIG A12, \
    ElMatrix_ ## SIG ABL, ElMatrix_ ## SIG ABR, \
    ElMatrix_ ## SIG A20, ElMatrix_ ## SIG A21, ElMatrix_ ## SIG A22 ); \
  ElError ElSlidePartitionUpDiagonalDist_ ## SIG \
  ( ElDistMatrix_ ## SIG ATL, ElDistMatrix_ ## SIG ATR, \
    ElDistMatrix_ ## SIG A00, ElDistMatrix_ ## SIG A01, ElDistMatrix_ ## SIG A02, \
    ElDistMatrix_ ## SIG A10, ElDistMatrix_ ## SIG A11, ElDistMatrix_ ## SIG A12, \
    ElDistMatrix_ ## SIG ABL, ElDistMatrix_ ## SIG ABR, \
    ElDistMatrix_ ## SIG A20, ElDistMatrix_ ## SIG A21, ElDistMatrix_ ## SIG A22 ); \
  ElError ElSlideLockedPartitionUpDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG ATL, ElMatrix_ ## SIG ATR, \
    ElConstMatrix_ ## SIG A00, ElConstMatrix_ ## SIG A01, ElConstMatrix_ ## SIG A02, \
    ElConstMatrix_ ## SIG A10, ElConstMatrix_ ## SIG A11, ElConstMatrix_ ## SIG A12, \
    ElMatrix_ ## SIG ABL, ElMatrix_ ## SIG ABR, \
    ElConstMatrix_ ## SIG A20, ElConstMatrix_ ## SIG A21, ElConstMatrix_ ## SIG A22 ); \
  ElError ElSlideLockedPartitionUpDiagonalDist_ ## SIG \
  ( ElDistMatrix_ ## SIG ATL, ElDistMatrix_ ## SIG ATR, \
    ElConstDistMatrix_ ## SIG A00, ElConstDistMatrix_ ## SIG A01, ElConstDistMatrix_ ## SIG A02, \
    ElConstDistMatrix_ ## SIG A10, ElConstDistMatrix_ ## SIG A11, ElConstDistMatrix_ ## SIG A12, \
    ElDistMatrix_ ## SIG ABL, ElDistMatrix_ ## SIG ABR, \
    ElConstDistMatrix_ ## SIG A20, ElConstDistMatrix_ ## SIG A21, ElConstDistMatrix_ ## SIG A22 );

#include "El/macros/CInstantiate.h"

} // extern "C"
