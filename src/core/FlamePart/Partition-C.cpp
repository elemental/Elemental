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
  /* Partition downwards from the top */ \
  ElError ElPartitionDown_ ## SIG \
  ( ElMatrix_ ## SIG A, \
    ElMatrix_ ## SIG AT, ElMatrix_ ## SIG AB, ElInt heightAT ) \
  { EL_TRY( PartitionDown( \
      *CReflect(A), *CReflect(AT), *CReflect(AB), heightAT ) ) } \
  ElError ElLockedPartitionDown_ ## SIG \
  ( ElConstMatrix_ ## SIG A, \
    ElMatrix_ ## SIG AT, ElMatrix_ ## SIG AB, ElInt heightAT ) \
  { EL_TRY( LockedPartitionDown( \
      *CReflect(A), *CReflect(AT), *CReflect(AB), heightAT ) ) } \
  ElError ElPartitionDownDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG AT, ElDistMatrix_ ## SIG AB, ElInt heightAT ) \
  { EL_TRY( PartitionDown( \
      *CReflect(A), *CReflect(AT), *CReflect(AB), heightAT ) ) } \
  ElError ElLockedPartitionDownDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG AT, ElDistMatrix_ ## SIG AB, ElInt heightAT ) \
  { EL_TRY( LockedPartitionDown( \
      *CReflect(A), *CReflect(AT), *CReflect(AB), heightAT ) ) } \
  /* Partition upwards from the bottom */ \
  ElError ElPartitionUp_ ## SIG \
  ( ElMatrix_ ## SIG A, \
    ElMatrix_ ## SIG AT, ElMatrix_ ## SIG AB, ElInt heightAB ) \
  { EL_TRY( PartitionUp( \
      *CReflect(A), *CReflect(AT), *CReflect(AB), heightAB ) ) } \
  ElError ElLockedPartitionUp_ ## SIG \
  ( ElConstMatrix_ ## SIG A, \
    ElMatrix_ ## SIG AT, ElMatrix_ ## SIG AB, ElInt heightAB ) \
  { EL_TRY( LockedPartitionUp( \
      *CReflect(A), *CReflect(AT), *CReflect(AB), heightAB ) ) } \
  ElError ElPartitionUpDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG AT, ElDistMatrix_ ## SIG AB, ElInt heightAB ) \
  { EL_TRY( PartitionUp( \
      *CReflect(A), *CReflect(AT), *CReflect(AB), heightAB ) ) } \
  ElError ElLockedPartitionUpDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG AT, ElDistMatrix_ ## SIG AB, ElInt heightAB ) \
  { EL_TRY( LockedPartitionUp( \
      *CReflect(A), *CReflect(AT), *CReflect(AB), heightAB ) ) } \
  /* Partition leftwards from the right */ \
  ElError ElPartitionLeft_ ## SIG \
  ( ElMatrix_ ## SIG A, \
    ElMatrix_ ## SIG AL, ElMatrix_ ## SIG AR, ElInt widthAR ) \
  { EL_TRY( PartitionLeft( \
      *CReflect(A), *CReflect(AL), *CReflect(AR), widthAR ) ) } \
  ElError ElPartitionLeftDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG AL, ElDistMatrix_ ## SIG AR, ElInt widthAR ) \
  { EL_TRY( PartitionLeft( \
      *CReflect(A), *CReflect(AL), *CReflect(AR), widthAR ) ) } \
  ElError ElLockedPartitionLeft_ ## SIG \
  ( ElConstMatrix_ ## SIG A, \
    ElMatrix_ ## SIG AL, ElMatrix_ ## SIG AR, ElInt widthAR ) \
  { EL_TRY( LockedPartitionLeft( \
      *CReflect(A), *CReflect(AL), *CReflect(AR), widthAR ) ) } \
  ElError ElLockedPartitionLeftDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG AL, ElDistMatrix_ ## SIG AR, ElInt widthAR ) \
  { EL_TRY( LockedPartitionLeft( \
      *CReflect(A), *CReflect(AL), *CReflect(AR), widthAR ) ) } \
  /* Partition rightwards from the left */ \
  ElError ElPartitionRight_ ## SIG \
  ( ElMatrix_ ## SIG A, \
    ElMatrix_ ## SIG AL, ElMatrix_ ## SIG AR, ElInt widthAL ) \
  { EL_TRY( PartitionRight( \
      *CReflect(A), *CReflect(AL), *CReflect(AR), widthAL ) ) } \
  ElError ElPartitionRightDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG AL, ElDistMatrix_ ## SIG AR, ElInt widthAL ) \
  { EL_TRY( PartitionRight( \
      *CReflect(A), *CReflect(AL), *CReflect(AR), widthAL ) ) } \
  ElError ElLockedPartitionRight_ ## SIG \
  ( ElConstMatrix_ ## SIG A, \
    ElMatrix_ ## SIG AL, ElMatrix_ ## SIG AR, ElInt widthAL ) \
  { EL_TRY( LockedPartitionRight( \
      *CReflect(A), *CReflect(AL), *CReflect(AR), widthAL ) ) } \
  ElError ElLockedPartitionRightDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG AL, ElDistMatrix_ ## SIG AR, ElInt widthAL ) \
  { EL_TRY( LockedPartitionRight( \
      *CReflect(A), *CReflect(AL), *CReflect(AR), widthAL ) ) } \
  /* Up the main diagonal */ \
  ElError ElPartitionUpDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG A, \
    ElMatrix_ ## SIG ATL, ElMatrix_ ## SIG ATR, \
    ElMatrix_ ## SIG ABL, ElMatrix_ ## SIG ABR, ElInt diagDist ) \
  { EL_TRY( PartitionUpDiagonal( \
      *CReflect(A), \
      *CReflect(ATL), *CReflect(ATR), *CReflect(ABL), *CReflect(ABR), \
      diagDist ) ) } \
  ElError ElPartitionUpDiagonalDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG ATL, ElDistMatrix_ ## SIG ATR, \
    ElDistMatrix_ ## SIG ABL, ElDistMatrix_ ## SIG ABR, ElInt diagDist ) \
  { EL_TRY( PartitionUpDiagonal( \
      *CReflect(A), \
      *CReflect(ATL), *CReflect(ATR), *CReflect(ABL), *CReflect(ABR), \
      diagDist ) ) } \
  ElError ElLockedPartitionUpDiagonal_ ## SIG \
  ( ElConstMatrix_ ## SIG A, \
    ElMatrix_ ## SIG ATL, ElMatrix_ ## SIG ATR, \
    ElMatrix_ ## SIG ABL, ElMatrix_ ## SIG ABR, ElInt diagDist ) \
  { EL_TRY( LockedPartitionUpDiagonal( \
      *CReflect(A), \
      *CReflect(ATL), *CReflect(ATR), *CReflect(ABL), *CReflect(ABR), \
      diagDist ) ) } \
  ElError ElLockedPartitionUpDiagonalDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG ATL, ElDistMatrix_ ## SIG ATR, \
    ElDistMatrix_ ## SIG ABL, ElDistMatrix_ ## SIG ABR, ElInt diagDist ) \
  { EL_TRY( LockedPartitionUpDiagonal( \
      *CReflect(A), \
      *CReflect(ATL), *CReflect(ATR), *CReflect(ABL), *CReflect(ABR), \
      diagDist ) ) } \
  /* Up an offset diagonal */ \
  ElError ElPartitionUpOffsetDiagonal_ ## SIG \
  ( ElInt offset, ElMatrix_ ## SIG A, \
    ElMatrix_ ## SIG ATL, ElMatrix_ ## SIG ATR, \
    ElMatrix_ ## SIG ABL, ElMatrix_ ## SIG ABR, ElInt diagDist ) \
  { EL_TRY( PartitionUpOffsetDiagonal( \
      offset, *CReflect(A), \
      *CReflect(ATL), *CReflect(ATR), *CReflect(ABL), *CReflect(ABR), \
      diagDist ) ) } \
  ElError ElPartitionUpOffsetDiagonalDist_ ## SIG \
  ( ElInt offset, ElDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG ATL, ElDistMatrix_ ## SIG ATR, \
    ElDistMatrix_ ## SIG ABL, ElDistMatrix_ ## SIG ABR, ElInt diagDist ) \
  { EL_TRY( PartitionUpOffsetDiagonal( \
      offset, *CReflect(A), \
      *CReflect(ATL), *CReflect(ATR), *CReflect(ABL), *CReflect(ABR), \
      diagDist ) ) } \
  ElError ElLockedPartitionUpOffsetDiagonal_ ## SIG \
  ( ElInt offset, ElConstMatrix_ ## SIG A, \
    ElMatrix_ ## SIG ATL, ElMatrix_ ## SIG ATR, \
    ElMatrix_ ## SIG ABL, ElMatrix_ ## SIG ABR, ElInt diagDist ) \
  { EL_TRY( LockedPartitionUpOffsetDiagonal( \
      offset, *CReflect(A), \
      *CReflect(ATL), *CReflect(ATR), *CReflect(ABL), *CReflect(ABR), \
      diagDist ) ) } \
  ElError ElLockedPartitionUpOffsetDiagonalDist_ ## SIG \
  ( ElInt offset, ElConstDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG ATL, ElDistMatrix_ ## SIG ATR, \
    ElDistMatrix_ ## SIG ABL, ElDistMatrix_ ## SIG ABR, ElInt diagDist ) \
  { EL_TRY( LockedPartitionUpOffsetDiagonal( \
      offset, *CReflect(A), \
      *CReflect(ATL), *CReflect(ATR), *CReflect(ABL), *CReflect(ABR), \
      diagDist ) ) } \
  /* Down the main diagonal */ \
  ElError ElPartitionDownDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG A, \
    ElMatrix_ ## SIG ATL, ElMatrix_ ## SIG ATR, \
    ElMatrix_ ## SIG ABL, ElMatrix_ ## SIG ABR, ElInt diagDist ) \
  { EL_TRY( PartitionDownDiagonal( \
      *CReflect(A), \
      *CReflect(ATL), *CReflect(ATR), *CReflect(ABL), *CReflect(ABR), \
      diagDist ) ) } \
  ElError ElPartitionDownDiagonalDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG ATL, ElDistMatrix_ ## SIG ATR, \
    ElDistMatrix_ ## SIG ABL, ElDistMatrix_ ## SIG ABR, ElInt diagDist ) \
  { EL_TRY( PartitionDownDiagonal( \
      *CReflect(A), \
      *CReflect(ATL), *CReflect(ATR), *CReflect(ABL), *CReflect(ABR), \
      diagDist ) ) } \
  ElError ElLockedPartitionDownDiagonal_ ## SIG \
  ( ElConstMatrix_ ## SIG A, \
    ElMatrix_ ## SIG ATL, ElMatrix_ ## SIG ATR, \
    ElMatrix_ ## SIG ABL, ElMatrix_ ## SIG ABR, ElInt diagDist ) \
  { EL_TRY( LockedPartitionDownDiagonal( \
      *CReflect(A), \
      *CReflect(ATL), *CReflect(ATR), *CReflect(ABL), *CReflect(ABR), \
      diagDist ) ) } \
  ElError ElLockedPartitionDownDiagonalDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG ATL, ElDistMatrix_ ## SIG ATR, \
    ElDistMatrix_ ## SIG ABL, ElDistMatrix_ ## SIG ABR, ElInt diagDist ) \
  { EL_TRY( LockedPartitionDownDiagonal( \
      *CReflect(A), \
      *CReflect(ATL), *CReflect(ATR), *CReflect(ABL), *CReflect(ABR), \
      diagDist ) ) } \
  /* Down an offset diagonal */ \
  ElError ElPartitionDownOffsetDiagonal_ ## SIG \
  ( ElInt offset, ElMatrix_ ## SIG A, \
    ElMatrix_ ## SIG ATL, ElMatrix_ ## SIG ATR, \
    ElMatrix_ ## SIG ABL, ElMatrix_ ## SIG ABR, ElInt diagDist ) \
  { EL_TRY( PartitionDownOffsetDiagonal( \
      offset, *CReflect(A), \
      *CReflect(ATL), *CReflect(ATR), *CReflect(ABL), *CReflect(ABR), \
      diagDist ) ) } \
  ElError ElPartitionDownOffsetDiagonalDist_ ## SIG \
  ( ElInt offset, ElDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG ATL, ElDistMatrix_ ## SIG ATR, \
    ElDistMatrix_ ## SIG ABL, ElDistMatrix_ ## SIG ABR, ElInt diagDist ) \
  { EL_TRY( PartitionDownOffsetDiagonal( \
      offset, *CReflect(A), \
      *CReflect(ATL), *CReflect(ATR), *CReflect(ABL), *CReflect(ABR), \
      diagDist ) ) } \
  ElError ElLockedPartitionDownOffsetDiagonal_ ## SIG \
  ( ElInt offset, ElConstMatrix_ ## SIG A, \
    ElMatrix_ ## SIG ATL, ElMatrix_ ## SIG ATR, \
    ElMatrix_ ## SIG ABL, ElMatrix_ ## SIG ABR, ElInt diagDist ) \
  { EL_TRY( LockedPartitionDownOffsetDiagonal( \
      offset, *CReflect(A), \
      *CReflect(ATL), *CReflect(ATR), *CReflect(ABL), *CReflect(ABR), \
      diagDist ) ) } \
  ElError ElLockedPartitionDownOffsetDiagonalDist_ ## SIG \
  ( ElInt offset, ElConstDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG ATL, ElDistMatrix_ ## SIG ATR, \
    ElDistMatrix_ ## SIG ABL, ElDistMatrix_ ## SIG ABR, ElInt diagDist ) \
  { EL_TRY( LockedPartitionDownOffsetDiagonal( \
      offset, *CReflect(A), \
      *CReflect(ATL), *CReflect(ATR), *CReflect(ABL), *CReflect(ABR), \
      diagDist ) ) }

#include "El/macros/CInstantiate.h"

} // extern "C"
