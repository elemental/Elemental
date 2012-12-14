/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

// To make our life easier. Undef'd at the bottom of the header
#define M  Matrix<T,Int>
#define DM DistMatrix<T,U,V,Int>

//
// PartitionUp
// 

template<typename T,typename Int>
void PartitionUp
( M& A, M& AT,
        M& AB, Int heightAB=Blocksize() );

template<typename T, Distribution U, Distribution V,typename Int>
void PartitionUp
( DM& A, DM& AT,
         DM& AB, Int heightAB=Blocksize() );

template<typename T,typename Int>
void LockedPartitionUp
( const M& A, M& AT,
              M& AB, Int heightAB=Blocksize() );

template<typename T, Distribution U, Distribution V,typename Int>
void LockedPartitionUp
( const DM& A, DM& AT,
               DM& AB, Int heightAB=Blocksize() );

//
// PartitionDown
//

template<typename T,typename Int>
void PartitionDown
( M& A, M& AT,
        M& AB, Int heightAT=Blocksize() );

template<typename T, Distribution U, Distribution V,typename Int>
void PartitionDown
( DM& A, DM& AT,
         DM& AB, Int heightAT=Blocksize() );

template<typename T,typename Int>
void LockedPartitionDown
( const M& A, M& AT,
              M& AB, Int heightAT=Blocksize() );

template<typename T, Distribution U, Distribution V,typename Int>
void LockedPartitionDown
( const DM& A, DM& AT,
               DM& AB, Int heightAT=Blocksize() );

//
// PartitionLeft
//

template<typename T,typename Int>
void PartitionLeft
( M& A, M& AL, M& AR, Int widthAR=Blocksize() );

template<typename T, Distribution U, Distribution V,typename Int>
void PartitionLeft
( DM& A, DM& AL, DM& AR, Int widthAR=Blocksize() );

template<typename T,typename Int>
void LockedPartitionLeft
( const M& A, M& AL, M& AR, Int widthAR=Blocksize() );

template<typename T, Distribution U, Distribution V,typename Int>
void LockedPartitionLeft
( const DM& A, DM& AL, DM& AR, Int widthAR=Blocksize() );

//
// PartitionRight
//

template<typename T,typename Int>
void PartitionRight
( M& A, M& AL, M& AR, Int widthAL=Blocksize() );

template<typename T, Distribution U, Distribution V,typename Int>
void PartitionRight
( DM& A, DM& AL, DM& AR, Int widthAL=Blocksize() );

template<typename T,typename Int>
void LockedPartitionRight
( const M& A, M& AL, M& AR, Int widthAL=Blocksize() );

template<typename T, Distribution U, Distribution V,typename Int>
void LockedPartitionRight
( const DM& A, DM& AL, DM& AR, Int widthAL=Blocksize() );

//
// PartitionUpDiagonal
//

template<typename T,typename Int>
void PartitionUpDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagABR=Blocksize() );

template<typename T, Distribution U, Distribution V,typename Int>
void PartitionUpDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagABR=Blocksize() );

template<typename T,typename Int>
void LockedPartitionUpDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagABR=Blocksize() );

template<typename T, Distribution U, Distribution V,typename Int>
void LockedPartitionUpDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagABR=Blocksize() );

//
// PartitionUpLeftDiagonal
//

template<typename T,typename Int>
void PartitionUpLeftDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagABR=Blocksize() );

template<typename T, Distribution U, Distribution V,typename Int>
void PartitionUpLeftDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagABR=Blocksize() );

template<typename T,typename Int>
void LockedPartitionUpLeftDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagABR=Blocksize() );

template<typename T, Distribution U, Distribution V,typename Int>
void LockedPartitionUpLeftDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagABR=Blocksize() );

//
// PartitionUpRightDiagonal
//

template<typename T,typename Int>
void PartitionUpRightDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagABR=Blocksize() );

template<typename T, Distribution U, Distribution V,typename Int>
void PartitionUpRightDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagABR=Blocksize() );

template<typename T,typename Int>
void LockedPartitionUpRightDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagABR=Blocksize() );

template<typename T, Distribution U, Distribution V,typename Int>
void LockedPartitionUpRightDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagABR=Blocksize() );

//
// PartitionDownDiagonal
//

template<typename T,typename Int>
void PartitionDownDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagATL=Blocksize() );

template<typename T, Distribution U, Distribution V,typename Int>
void PartitionDownDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagATL=Blocksize() );

template<typename T,typename Int>
void LockedPartitionDownDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagATL=Blocksize() );

template<typename T, Distribution U, Distribution V,typename Int>
void LockedPartitionDownDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagATL=Blocksize() );

//
// PartitionDownLeftDiagonal
//

template<typename T,typename Int>
void PartitionDownLeftDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagATL=Blocksize() );

template<typename T, Distribution U, Distribution V,typename Int>
void PartitionDownLeftDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagATL=Blocksize() );

template<typename T,typename Int>
void LockedPartitionDownLeftDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagATL=Blocksize() );

template<typename T, Distribution U, Distribution V,typename Int>
void LockedPartitionDownLeftDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagATL=Blocksize() );

//
// PartitionDownRightDiagonal
//

template<typename T,typename Int>
void PartitionDownRightDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagATL=Blocksize() );

template<typename T, Distribution U, Distribution V,typename Int>
void PartitionDownRightDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagATL=Blocksize() );

template<typename T,typename Int>
void LockedPartitionDownRightDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagATL=Blocksize() );

template<typename T, Distribution U, Distribution V,typename Int>
void LockedPartitionDownRightDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagATL=Blocksize() );

#undef DM
#undef M

} // namespace elem

