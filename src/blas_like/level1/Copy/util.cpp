/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace copy {
namespace util {

template<typename T>
void InterleaveMatrix
( Int height, Int width,
  const T* A, Int colStrideA, Int rowStrideA,
        T* B, Int colStrideB, Int rowStrideB )
{
    // TODO: Add OpenMP parallelization and/or optimize
    if( colStrideA == 1 && colStrideB == 1 )
    {
        for( Int j=0; j<width; ++j )
            MemCopy( &B[j*rowStrideB], &A[j*rowStrideA], height );
    }
    else
    {
        for( Int j=0; j<width; ++j )
            StridedMemCopy
            ( &B[j*rowStrideB], colStrideB,
              &A[j*rowStrideA], colStrideA, height );
    }
}

// TODO: ColStridedPack

template<typename T>
void ColStridedUnpack
( Int height, Int width,
  Int colAlign, Int colStride,
  const T* APortions, Int portionSize,
        T* B,         Int BLDim )
{
    for( Int k=0; k<colStride; ++k )
    {
        const Int colShift = Shift_( k, colAlign, colStride );
        const Int localHeight = Length_( height, colShift, colStride );
        InterleaveMatrix
        ( localHeight, width,
          &APortions[k*portionSize], 1,         localHeight,
          &B[colShift],              colStride, BLDim );
    }
}

template<typename T>
void PartialColStridedUnpack
( Int height, Int width,
  Int colAlign, Int colStride,
  Int colStrideUnion, Int colStridePart, Int colRankPart,
  Int colShiftB,
  const T* APortions, Int portionSize,
        T* B,         Int BLDim )
{
    for( Int k=0; k<colStrideUnion; ++k )
    {
        const Int colShift =
            Shift_( colRankPart+k*colStridePart, colAlign, colStride );
        const Int colOffset = (colShift-colShiftB) / colStridePart;
        const Int localHeight = Length_( height, colShift, colStride );
        InterleaveMatrix
        ( localHeight, width,
          &APortions[k*portionSize], 1,              localHeight,
          &B[colOffset],             colStrideUnion, BLDim );
    }
}

template<typename T>
void RowStridedPack
( Int height, Int width,
  Int rowAlign, Int rowStride,
  const T* A,         Int ALDim,
        T* BPortions, Int portionSize )
{
    for( Int k=0; k<rowStride; ++k )
    {
        const Int rowShift = Shift_( k, rowAlign, rowStride );
        const Int localWidth = Length_( width, rowShift, rowStride );
        InterleaveMatrix
        ( height, localWidth,
          &A[rowShift*ALDim],        1, rowStride*ALDim,
          &BPortions[k*portionSize], 1, height );
    }
}

template<typename T>
void RowStridedUnpack
( Int height, Int width,
  Int rowAlign, Int rowStride,
  const T* APortions, Int portionSize,
        T* B,         Int BLDim )
{
    for( Int k=0; k<rowStride; ++k )
    {
        const Int rowShift = Shift_( k, rowAlign, rowStride );
        const Int localWidth = Length_( width, rowShift, rowStride );
        InterleaveMatrix
        ( height, localWidth,
          &APortions[k*portionSize], 1, height,
          &B[rowShift*BLDim],        1, rowStride*BLDim );
    }
}

// TODO: PartialRowStridedPack

template<typename T>
void PartialRowStridedUnpack
( Int height, Int width,
  Int rowAlign, Int rowStride,
  Int rowStrideUnion, Int rowStridePart, Int rowRankPart,
  Int rowShiftB,
  const T* APortions, Int portionSize,
        T* B,         Int BLDim )
{
    for( Int k=0; k<rowStrideUnion; ++k )
    {
        const Int rowShift =
            Shift_( rowRankPart+k*rowStridePart, rowAlign, rowStride );
        const Int rowOffset = (rowShift-rowShiftB) / rowStridePart;
        const Int localWidth = Length_( width, rowShift, rowStride );
        InterleaveMatrix
        ( height, localWidth,
          &APortions[k*portionSize], 1, height,
          &B[rowOffset*BLDim],       1, rowStrideUnion*BLDim );
    }
}

// NOTE: This is implicitly column-major
template<typename T>
void StridedUnpack
( Int height, Int width,
  Int colAlign, Int colStride,
  Int rowAlign, Int rowStride,
  const T* APortions, Int portionSize,
        T* B,         Int BLDim )
{
    for( Int l=0; l<rowStride; ++l )
    {
        const Int rowShift = Shift_( l, rowAlign, rowStride );
        const Int localWidth = Length_( width, rowShift, rowStride );
        for( Int k=0; k<colStride; ++k )
        {
            const Int colShift = Shift_( k, colAlign, colStride );
            const Int localHeight = Length_( height, colShift, colStride );
            InterleaveMatrix
            ( localHeight, localWidth,
              &APortions[(k+l*colStride)*portionSize], 1, localHeight,
              &B[colShift+rowShift*BLDim], colStride, rowStride*BLDim );
        }
    }
}

#define PROTO(T) \
  template void InterleaveMatrix \
  ( Int height, Int width, \
    const T* A, Int colStrideA, Int rowStrideA, \
          T* B, Int colStrideB, Int rowStrideB ); \
  template void ColStridedUnpack \
  ( Int height, Int width, \
    Int colAlign, Int colStride, \
    const T* APortions, Int portionSize, \
          T* B,         Int BLDim ); \
  template void PartialColStridedUnpack \
  ( Int height, Int width, \
    Int colAlign, Int colStride, \
    Int colStrideUnion, Int colStridePart, Int colRankPart, \
    Int colShiftB, \
    const T* APortions, Int portionSize, \
          T* B,         Int BLDim ); \
  template void RowStridedPack \
  ( Int height, Int width, \
    Int rowAlign, Int rowStride, \
    const T* A,         Int ALDim, \
          T* BPortions, Int portionSize ); \
  template void RowStridedUnpack \
  ( Int height, Int width, \
    Int rowAlign, Int rowStride, \
    const T* APortions, Int portionSize, \
          T* B,         Int BLDim ); \
  template void PartialRowStridedUnpack \
  ( Int height, Int width, \
    Int rowAlign, Int rowStride, \
    Int rowStrideUnion, Int rowStridePart, Int rowRankPart, \
    Int rowShiftB, \
    const T* APortions, Int portionSize, \
          T* B,         Int BLDim ); \
  template void StridedUnpack \
  ( Int height, Int width, \
    Int colAlign, Int colStride, \
    Int rowAlign, Int rowStride, \
    const T* APortions, Int portionSize, \
          T* B,         Int BLDim );

#include "El/macros/Instantiate.h"

} // namespace util
} // namespace copy
} // namespace El
