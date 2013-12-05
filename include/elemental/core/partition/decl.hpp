/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_PARTITION_DECL_HPP
#define ELEM_CORE_PARTITION_DECL_HPP

namespace elem {

// To make our life easier. Undef'd at the bottom of the header
#define M  Matrix<T>
#define DM DistMatrix<T,U,V>

//
// PartitionUp
// 

template<typename T>
void PartitionUp
( M& A, M& AT,
        M& AB, Int heightAB=Blocksize() );

template<typename T, Distribution U, Distribution V>
void PartitionUp
( DM& A, DM& AT,
         DM& AB, Int heightAB=Blocksize() );

template<typename T>
void LockedPartitionUp
( const M& A, M& AT,
              M& AB, Int heightAB=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedPartitionUp
( const DM& A, DM& AT,
               DM& AB, Int heightAB=Blocksize() );

//
// PartitionDown
//

template<typename T>
void PartitionDown
( M& A, M& AT,
        M& AB, Int heightAT=Blocksize() );

template<typename T, Distribution U, Distribution V>
void PartitionDown
( DM& A, DM& AT,
         DM& AB, Int heightAT=Blocksize() );

template<typename T>
void LockedPartitionDown
( const M& A, M& AT,
              M& AB, Int heightAT=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedPartitionDown
( const DM& A, DM& AT,
               DM& AB, Int heightAT=Blocksize() );

//
// PartitionLeft
//

template<typename T>
void PartitionLeft
( M& A, M& AL, M& AR, Int widthAR=Blocksize() );

template<typename T, Distribution U, Distribution V>
void PartitionLeft
( DM& A, DM& AL, DM& AR, Int widthAR=Blocksize() );

template<typename T>
void LockedPartitionLeft
( const M& A, M& AL, M& AR, Int widthAR=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedPartitionLeft
( const DM& A, DM& AL, DM& AR, Int widthAR=Blocksize() );

//
// PartitionRight
//

template<typename T>
void PartitionRight
( M& A, M& AL, M& AR, Int widthAL=Blocksize() );

template<typename T, Distribution U, Distribution V>
void PartitionRight
( DM& A, DM& AL, DM& AR, Int widthAL=Blocksize() );

template<typename T>
void LockedPartitionRight
( const M& A, M& AL, M& AR, Int widthAL=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedPartitionRight
( const DM& A, DM& AL, DM& AR, Int widthAL=Blocksize() );

//
// PartitionUpDiagonal
//

template<typename T>
void PartitionUpDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagDist=Blocksize() );

template<typename T, Distribution U, Distribution V>
void PartitionUpDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagDist=Blocksize() );

template<typename T>
void LockedPartitionUpDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagDist=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedPartitionUpDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagDist=Blocksize() );

//
// PartitionUpOffsetDiagonal
//

template<typename T>
void PartitionUpOffsetDiagonal
( Int offset,
  M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagDist=Blocksize() );

template<typename T, Distribution U, Distribution V>
void PartitionUpOffsetDiagonal
( Int offset,
  DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagDist=Blocksize() );

template<typename T>
void LockedPartitionUpOffsetDiagonal
( Int offset,
  const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagDist=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedPartitionUpOffsetDiagonal
( Int offset,
  const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagDist=Blocksize() );

//
// PartitionDownDiagonal
//

template<typename T>
void PartitionDownDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagDist=Blocksize() );

template<typename T, Distribution U, Distribution V>
void PartitionDownDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagDist=Blocksize() );

template<typename T>
void LockedPartitionDownDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagDist=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedPartitionDownDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagDist=Blocksize() );

//
// PartitionDownOffsetDiagonal
//

template<typename T>
void PartitionDownOffsetDiagonal
( Int offset,
  M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagDist=Blocksize() );

template<typename T, Distribution U, Distribution V>
void PartitionDownOffsetDiagonal
( Int offset,
  DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagDist=Blocksize() );

template<typename T>
void LockedPartitionDownOffsetDiagonal
( Int offset,
  const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagDist=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedPartitionDownOffsetDiagonal
( Int offset,
  const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagDist=Blocksize() );

#undef DM
#undef M

} // namespace elem

#endif // ifndef ELEM_CORE_PARTITION_DECL_HPP
