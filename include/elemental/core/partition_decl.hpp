/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
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

