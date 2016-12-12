/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_COPY_UTIL_HPP
#define EL_BLAS_COPY_UTIL_HPP

namespace El {
namespace copy {
namespace util {

template<typename T>
void InterleaveMatrix
( Int height, Int width,
  const T* A, Int colStrideA, Int rowStrideA,
        T* B, Int colStrideB, Int rowStrideB )
{
    if( colStrideA == 1 && colStrideB == 1 )
    {
        lapack::Copy( 'F', height, width, A, rowStrideA, B, rowStrideB );
    }
    else
    {
#ifdef EL_HAVE_MKL
        mkl::omatcopy
        ( NORMAL, height, width, T(1),
          A, rowStrideA, colStrideA,
          B, rowStrideB, colStrideB );
#else
        for( Int j=0; j<width; ++j )
            StridedMemCopy
            ( &B[j*rowStrideB], colStrideB,
              &A[j*rowStrideA], colStrideA, height );
#endif
    }
}

template<typename T>
void ColStridedPack
( Int height, Int width,
  Int colAlign, Int colStride,
  const T* A,         Int ALDim,
        T* BPortions, Int portionSize )
{
    for( Int k=0; k<colStride; ++k )
    {
        const Int colShift = Shift_( k, colAlign, colStride );
        const Int localHeight = Length_( height, colShift, colStride );
        InterleaveMatrix
        ( localHeight, width,
          &A[colShift],              colStride, ALDim,
          &BPortions[k*portionSize], 1,         localHeight );
    }
}

// TODO(poulson): Use this routine
template<typename T>
void ColStridedColumnPack
( Int height,
  Int colAlign, Int colStride,
  const T* A,
        T* BPortions, Int portionSize )
{
    for( Int k=0; k<colStride; ++k )
    {
        const Int colShift = Shift_( k, colAlign, colStride );
        const Int localHeight = Length_( height, colShift, colStride );
        StridedMemCopy
        ( &BPortions[k*portionSize], 1,
          &A[colShift],              colStride, localHeight );
    }
}

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
void BlockedColStridedUnpack
( Int height, Int width,
  Int colAlign, Int colStride,
  Int blockHeight, Int colCut,
  const T* APortions, Int portionSize,
        T* B,         Int BLDim )
{
    const Int firstBlockHeight = blockHeight - colCut;
    for( Int portion=0; portion<colStride; ++portion )
    {
        const T* APortion = &APortions[portion*portionSize];
        const Int colShift = Shift_( portion, colAlign, colStride );
        const Int localHeight =
          BlockedLength_( height, colShift, blockHeight, colCut, colStride );

        // Loop over the block rows from this portion
        Int blockRow = colShift;
        Int rowIndex =
          ( colShift==0 ? 0 : firstBlockHeight + (colShift-1)*blockHeight );
        Int packedRowIndex = 0;
        while( rowIndex < height )
        {
            const Int thisBlockHeight =
              ( blockRow == 0 ?
                firstBlockHeight :
                Min(blockHeight,height-rowIndex) );

            lapack::Copy
            ( 'F', thisBlockHeight, width,
              &APortion[packedRowIndex], localHeight,
              &B[rowIndex],              BLDim );

            blockRow += colStride;
            rowIndex += thisBlockHeight + (colStride-1)*blockHeight;
            packedRowIndex += thisBlockHeight;
        }
    }
}

template<typename T>
void PartialColStridedPack
( Int height, Int width,
  Int colAlign, Int colStride,
  Int colStrideUnion, Int colStridePart, Int colRankPart,
  Int colShiftA,
  const T* A,         Int ALDim,
        T* BPortions, Int portionSize )
{
    for( Int k=0; k<colStrideUnion; ++k )
    {
        const Int colShift =
            Shift_( colRankPart+k*colStridePart, colAlign, colStride );
        const Int colOffset = (colShift-colShiftA) / colStridePart;
        const Int localHeight = Length_( height, colShift, colStride );
        InterleaveMatrix
        ( localHeight, width,
          &A[colOffset],             colStrideUnion, ALDim,
          &BPortions[k*portionSize], 1,              localHeight );
    }
}

template<typename T>
void PartialColStridedColumnPack
( Int height,
  Int colAlign, Int colStride,
  Int colStrideUnion, Int colStridePart, Int colRankPart,
  Int colShiftA,
  const T* A,
        T* BPortions, Int portionSize )
{
    for( Int k=0; k<colStrideUnion; ++k )
    {
        const Int colShift =
            Shift_( colRankPart+k*colStridePart, colAlign, colStride );
        const Int colOffset = (colShift-colShiftA) / colStridePart;
        const Int localHeight = Length_( height, colShift, colStride );
        StridedMemCopy
        ( &BPortions[k*portionSize], 1,
          &A[colOffset],             colStrideUnion, localHeight );
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
void PartialColStridedColumnUnpack
( Int height,
  Int colAlign, Int colStride,
  Int colStrideUnion, Int colStridePart, Int colRankPart,
  Int colShiftB,
  const T* APortions, Int portionSize,
        T* B )
{
    for( Int k=0; k<colStrideUnion; ++k )
    {
        const Int colShift =
            Shift_( colRankPart+k*colStridePart, colAlign, colStride );
        const Int colOffset = (colShift-colShiftB) / colStridePart;
        const Int localHeight = Length_( height, colShift, colStride );
        StridedMemCopy
        ( &B[colOffset],             colStrideUnion,
          &APortions[k*portionSize], 1,              localHeight );
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
        lapack::Copy
        ( 'F', height, localWidth,
          &A[rowShift*ALDim],        rowStride*ALDim,
          &BPortions[k*portionSize], height );
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
        lapack::Copy
        ( 'F', height, localWidth,
          &APortions[k*portionSize], height,
          &B[rowShift*BLDim],        rowStride*BLDim );
    }
}

template<typename T>
void BlockedRowStridedUnpack
( Int height, Int width,
  Int rowAlign, Int rowStride,
  Int blockWidth, Int rowCut,
  const T* APortions, Int portionSize,
        T* B,         Int BLDim )
{
    const Int firstBlockWidth = blockWidth - rowCut;
    for( Int portion=0; portion<rowStride; ++portion )
    {
        const T* APortion = &APortions[portion*portionSize];
        const Int rowShift = Shift_( portion, rowAlign, rowStride );
        // Loop over the block columns from this portion
        Int blockCol = rowShift;
        Int colIndex =
          ( rowShift==0 ? 0 : firstBlockWidth + (rowShift-1)*blockWidth );
        Int packedColIndex = 0;

        while( colIndex < width )
        {
            const Int thisBlockWidth =
              ( blockCol == 0 ?
                firstBlockWidth :
                Min(blockWidth,width-colIndex) );

            lapack::Copy
            ( 'F', height, thisBlockWidth,
              &APortion[packedColIndex*height], height,
              &B[colIndex*BLDim],               BLDim );

            blockCol += rowStride;
            colIndex += thisBlockWidth + (rowStride-1)*blockWidth;
            packedColIndex += thisBlockWidth;
        }
    }
}

template<typename T>
void BlockedRowFilter
( Int height, Int width,
  Int rowShift, Int rowStride,
  Int blockWidth, Int rowCut,
  const T* A, Int ALDim,
        T* B, Int BLDim )
{
    EL_DEBUG_CSE
    const Int firstBlockWidth = blockWidth - rowCut;

    // Loop over the block columns from this portion
    Int blockCol = rowShift;
    Int colIndex =
      ( rowShift==0 ? 0 : firstBlockWidth + (rowShift-1)*blockWidth );
    Int packedColIndex = 0;

    while( colIndex < width )
    {
        const Int thisBlockWidth =
          ( blockCol == 0 ?
            firstBlockWidth :
            Min(blockWidth,width-colIndex) );

        lapack::Copy
        ( 'F', height, thisBlockWidth,
          &A[colIndex      *ALDim], ALDim,
          &B[packedColIndex*BLDim], BLDim );

        blockCol += rowStride;
        colIndex += thisBlockWidth + (rowStride-1)*blockWidth;
        packedColIndex += thisBlockWidth;
    }
}

template<typename T>
void BlockedColFilter
( Int height, Int width,
  Int colShift, Int colStride,
  Int blockHeight, Int colCut,
  const T* A, Int ALDim,
        T* B, Int BLDim )
{
    EL_DEBUG_CSE
    const Int firstBlockHeight = blockHeight - colCut;

    // Loop over the block rows from this portion
    Int blockRow = colShift;
    Int rowIndex =
      ( colShift==0 ? 0 : firstBlockHeight + (colShift-1)*blockHeight );
    Int packedRowIndex = 0;

    while( rowIndex < height )
    {
        const Int thisBlockHeight =
          ( blockRow == 0 ?
            firstBlockHeight :
            Min(blockHeight,height-rowIndex) );

        lapack::Copy
        ( 'F', thisBlockHeight, width,
          &A[rowIndex],       ALDim,
          &B[packedRowIndex], BLDim );

        blockRow += colStride;
        rowIndex += thisBlockHeight + (colStride-1)*blockHeight;
        packedRowIndex += thisBlockHeight;
    }
}

template<typename T>
void PartialRowStridedPack
( Int height, Int width,
  Int rowAlign, Int rowStride,
  Int rowStrideUnion, Int rowStridePart, Int rowRankPart,
  Int rowShiftA,
  const T* A,         Int ALDim,
        T* BPortions, Int portionSize )
{
    for( Int k=0; k<rowStrideUnion; ++k )
    {
        const Int rowShift =
            Shift_( rowRankPart+k*rowStridePart, rowAlign, rowStride );
        const Int rowOffset = (rowShift-rowShiftA) / rowStridePart;
        const Int localWidth = Length_( width, rowShift, rowStride );
        lapack::Copy
        ( 'F', height, localWidth,
          &A[rowOffset*ALDim],       rowStrideUnion*ALDim,
          &BPortions[k*portionSize], height );
    }
}
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
        lapack::Copy
        ( 'F', height, localWidth,
          &APortions[k*portionSize], height,
          &B[rowOffset*BLDim],       rowStrideUnion*BLDim );
    }
}

// NOTE: This is implicitly column-major
template<typename T>
void StridedPack
( Int height, Int width,
  Int colAlign, Int colStride,
  Int rowAlign, Int rowStride,
  const T* A,         Int ALDim,
        T* BPortions, Int portionSize )
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
              &A[colShift+rowShift*ALDim], colStride, rowStride*ALDim,
              &BPortions[(k+l*colStride)*portionSize], 1, localHeight );
        }
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

} // namespace util
} // namespace copy
} // namespace El

#endif // ifndef EL_BLAS_COPY_UTIL_HPP
