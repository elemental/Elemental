/*
   Copyright (c) 2009-2011, Jack Poulson
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
#ifndef ELEMENTAL_PARTITIONING_HPP
#define ELEMENTAL_PARTITIONING_HPP 1

#include "elemental/core/dist_matrix.hpp"

// To make our life easier. Undef'd at the bottom of the header
#define M  Matrix<T>
#define DM DistMatrix<T,U,V>

namespace elemental {

template<typename T>
void PartitionUp
( M& A, M& AT,
        M& AB, int heightAB=Blocksize() );

template<typename T, Distribution U, Distribution V>
void PartitionUp
( DM& A, DM& AT,
         DM& AB, int heightAB=Blocksize() );

template<typename T>
void PartitionDown
( M& A, M& AT,
        M& AB, int heightAT=Blocksize() );

template<typename T, Distribution U, Distribution V>
void PartitionDown
( DM& A, DM& AT,
         DM& AB, int heightAT=Blocksize() );

template<typename T>
void PartitionLeft
( M& A, M& AL, M& AR, int widthAR=Blocksize() );

template<typename T, Distribution U, Distribution V>
void PartitionLeft
( DM& A, DM& AL, DM& AR, int widthAR=Blocksize() );

template<typename T>
void PartitionRight
( M& A, M& AL, M& AR, int widthAL=Blocksize() );

template<typename T, Distribution U, Distribution V>
void PartitionRight
( DM& A, DM& AL, DM& AR, int widthAL=Blocksize() );

template<typename T>
void PartitionUpDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, int diagABR=Blocksize() );

template<typename T, Distribution U, Distribution V>
void PartitionUpDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, int diagABR=Blocksize() );

template<typename T>
void PartitionUpLeftDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, int diagABR=Blocksize() );

template<typename T, Distribution U, Distribution V>
void PartitionUpLeftDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, int diagABR=Blocksize() );

template<typename T>
void PartitionUpRightDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, int diagABR=Blocksize() );

template<typename T, Distribution U, Distribution V>
void PartitionUpRightDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, int diagABR=Blocksize() );

template<typename T>
void PartitionDownDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, int diagATL=Blocksize() );

template<typename T, Distribution U, Distribution V>
void PartitionDownDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, int diagATL=Blocksize() );

template<typename T>
void PartitionDownLeftDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, int diagATL=Blocksize() );

template<typename T, Distribution U, Distribution V>
void PartitionDownLeftDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, int diagATL=Blocksize() );

template<typename T>
void PartitionDownRightDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, int diagATL=Blocksize() );

template<typename T, Distribution U, Distribution V>
void PartitionDownRightDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, int diagATL=Blocksize() );

template<typename T>
void LockedPartitionUp
( const M& A, M& AT,
              M& AB, int heightAB=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedPartitionUp
( const DM& A, DM& AT,
               DM& AB, int heightAB=Blocksize() );

template<typename T>
void LockedPartitionDown
( const M& A, M& AT,
              M& AB, int heightAT=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedPartitionDown
( const DM& A, DM& AT,
               DM& AB, int heightAT=Blocksize() );

template<typename T>
void LockedPartitionLeft
( const M& A, M& AL, M& AR, int widthAR=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedPartitionLeft
( const DM& A, DM& AL, DM& AR, int widthAR=Blocksize() );

template<typename T>
void LockedPartitionRight
( const M& A, M& AL, M& AR, int widthAL=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedPartitionRight
( const DM& A, DM& AL, DM& AR, int widthAL=Blocksize() );

template<typename T>
void LockedPartitionUpDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, int diagABR=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedPartitionUpDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, int diagABR=Blocksize() );

template<typename T>
void LockedPartitionUpLeftDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, int diagABR=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedPartitionUpLeftDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, int diagABR=Blocksize() );

template<typename T>
void LockedPartitionUpRightDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, int diagABR=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedPartitionUpRightDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, int diagABR=Blocksize() );

template<typename T>
void LockedPartitionDownDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, int diagATL=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedPartitionDownDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, int diagATL=Blocksize() );

template<typename T>
void LockedPartitionDownLeftDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, int diagATL=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedPartitionDownLeftDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, int diagATL=Blocksize() );

template<typename T>
void LockedPartitionDownRightDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, int diagATL=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedPartitionDownRightDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, int diagATL=Blocksize() );

template<typename T>
void RepartitionUp
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2, int bsize=Blocksize() );

template<typename T, Distribution U, Distribution V>
void RepartitionUp
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2, int bsize=Blocksize() );

template<typename T>
void RepartitionDown
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2, int bsize=Blocksize() );

template<typename T, Distribution U, Distribution V>
void RepartitionDown
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2, int bsize=Blocksize() );

template<typename T>
void RepartitionLeft
( M& AL, M& AR,
  M& A0, M& A1, M& A2, int bsize=Blocksize() );

template<typename T, Distribution U, Distribution V>
void RepartitionLeft
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2, int bsize=Blocksize() );

template<typename T>
void RepartitionRight
( M& AL, M& AR,
  M& A0, M& A1, M& A2, int bsize=Blocksize() );

template<typename T, Distribution U, Distribution V>
void RepartitionRight
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2, int bsize=Blocksize() );

template<typename T>
void RepartitionUpDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22, int bsize=Blocksize() );

template<typename T, Distribution U, Distribution V>
void RepartitionUpDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22, int bsize=Blocksize() );

template<typename T>
void RepartitionDownDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22, int bsize=Blocksize() );

template<typename T, Distribution U, Distribution V>
void RepartitionDownDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22, int bsize=Blocksize() );

template<typename T>
void LockedRepartitionUp
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2, int bsize=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedRepartitionUp
( const DM& AT, DM& A0,
                DM& A1,
  const DM& AB, DM& A2, int bsize=Blocksize() );

template<typename T>
void LockedRepartitionDown
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2, int bsize=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedRepartitionDown
( const DM& AT, DM& A0,
                DM& A1,
  const DM& AB, DM& A2, int bsize=Blocksize() );

template<typename T>
void LockedRepartitionLeft
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2, int bsize=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedRepartitionLeft
( const DM& AL, const DM& AR,
  DM& A0, DM& A1, DM& A2, int bsize=Blocksize() );

template<typename T>
void LockedRepartitionRight
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2, int bsize=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedRepartitionRight
( const DM& AL, const DM& AR,
  DM& A0, DM& A1, DM& A2, int bsize=Blocksize() );

template<typename T>
void LockedRepartitionUpDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22, int bsize=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedRepartitionUpDiagonal
( const DM& ATL, const DM& ATR, DM& A00, DM& A01, DM& A02,
                                DM& A10, DM& A11, DM& A12,
  const DM& ABL, const DM& ABR, DM& A20, DM& A21, DM& A22, 
  int bsize=Blocksize() );

template<typename T>
void LockedRepartitionDownDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22, int bsize=Blocksize() );

template<typename T, Distribution U, Distribution V>
void LockedRepartitionDownDiagonal
( const DM& ATL, const DM& ATR, DM& A00, DM& A01, DM& A02,
                                DM& A10, DM& A11, DM& A12,
  const DM& ABL, const DM& ABR, DM& A20, DM& A21, DM& A22, 
  int bsize=Blocksize() );

template<typename T>
void SlidePartitionUp
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2 );

template<typename T, Distribution U, Distribution V>
void SlidePartitionUp
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2 );

template<typename T>
void SlidePartitionDown
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2 );

template<typename T, Distribution U, Distribution V>
void SlidePartitionDown
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2 );

template<typename T>
void SlidePartitionLeft
( M& AL, M& AR,
  M& A0, M& A1, M& A2 );

template<typename T, Distribution U, Distribution V>
void SlidePartitionLeft
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2 );

template<typename T>
void SlidePartitionRight
( M& AL, M& AR,
  M& A0, M& A1, M& A2 );

template<typename T, Distribution U, Distribution V>
void SlidePartitionRight
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2 );

template<typename T>
void SlidePartitionUpDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22 );

template<typename T, Distribution U, Distribution V>
void SlidePartitionUpDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22 );

template<typename T>
void SlidePartitionDownDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22 );

template<typename T, Distribution U, Distribution V>
void SlidePartitionDownDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22 );

template<typename T>
void SlideLockedPartitionUp
( M& AT, const M& A0,
         const M& A1,
  M& AB, const M& A2 );

template<typename T, Distribution U, Distribution V>
void SlideLockedPartitionUp
( DM& AT, const DM& A0,
          const DM& A1,
  DM& AB, const DM& A2 );

template<typename T>
void SlideLockedPartitionDown
( M& AT, const M& A0,
         const M& A1,
  M& AB, const M& A2 );

template<typename T, Distribution U, Distribution V>
void SlideLockedPartitionDown
( DM& AT, const DM& A0,
          const DM& A1,
  DM& AB, const DM& A2 );

template<typename T>
void SlideLockedPartitionLeft
( M& AL, M& AR,
  const M& A0, const M& A1, const M& A2 ); 

template<typename T, Distribution U, Distribution V>
void SlideLockedPartitionLeft
( DM& AL, DM& AR,
  const DM& A0, const DM& A1, const DM& A2 );

template<typename T>
void SlideLockedPartitionRight
( M& AL, M& AR,
  const M& A0, const M& A1, const M& A2 );

template<typename T, Distribution U, Distribution V>
void SlideLockedPartitionRight
( DM& AL, DM& AR,
  const DM& A0, const DM& A1, const DM& A2 );

template<typename T>
void SlideLockedPartitionUpDiagonal
( M& ATL, M& ATR, const M& A00, const M& A01, const M& A02,
                  const M& A10, const M& A11, const M& A12,
  M& ABL, M& ABR, const M& A20, const M& A21, const M& A22 );

template<typename T, Distribution U, Distribution V>
void SlideLockedPartitionUpDiagonal
( DM& ATL, DM& ATR, const DM& A00, const DM& A01, const DM& A02,
                    const DM& A10, const DM& A11, const DM& A12,
  DM& ABL, DM& ABR, const DM& A20, const DM& A21, const DM& A22 );

template<typename T>
void SlideLockedPartitionDownDiagonal
( M& ATL, M& ATR, const M& A00, const M& A01, const M& A02,
                  const M& A10, const M& A11, const M& A12,
  M& ABL, M& ABR, const M& A20, const M& A21, const M& A22 );

template<typename T, Distribution U, Distribution V>
void SlideLockedPartitionDownDiagonal
( DM& ATL, DM& ATR, const DM& A00, const DM& A01, const DM& A02,
                    const DM& A10, const DM& A11, const DM& A12,
  DM& ABL, DM& ABR, const DM& A20, const DM& A21, const DM& A22 );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename T>
inline void
PartitionUp
( M& A, M& AT,
        M& AB, int heightAB )
{
#ifndef RELEASE
    PushCallStack("PartitionUp [Matrix]");
    if( heightAB < 0 )
        throw std::logic_error
        ("Height of bottom partition must be non-negative");
#endif
    heightAB = std::min(heightAB,A.Height());
    const int heightAT = A.Height()-heightAB;
    AT.View( A, 0,        0, heightAT, A.Width() );
    AB.View( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
PartitionUp
( DM& A, DM& AT,
         DM& AB, int heightAB )
{
#ifndef RELEASE
    PushCallStack("PartitionUp [DistMatrix]");
    if( heightAB < 0 )
        throw std::logic_error
        ("Height of bottom partition must be non-negative");
#endif
    heightAB = std::min(heightAB,A.Height());
    const int heightAT = A.Height()-heightAB;
    AT.View( A, 0,        0, heightAT, A.Width() );
    AB.View( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
PartitionDown
( M& A, M& AT,
        M& AB, int heightAT ) 
{
#ifndef RELEASE
    PushCallStack("PartitionDown [Matrix]");
    if( heightAT < 0 )
        throw std::logic_error("Height of top partition must be non-negative");
#endif
    heightAT = std::min(heightAT,A.Height());
    const int heightAB = A.Height()-heightAT;
    AT.View( A, 0,        0, heightAT, A.Width() );
    AB.View( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
PartitionDown
( DM& A, DM& AT,
         DM& AB, int heightAT )
{
#ifndef RELEASE
    PushCallStack("PartitionDown [DistMatrix]");
    if( heightAT < 0 )
        throw std::logic_error("Height of top partition must be non-negative");
#endif
    heightAT = std::min(heightAT,A.Height());
    const int heightAB = A.Height()-heightAT;
    AT.View( A, 0,        0, heightAT, A.Width() );
    AB.View( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
PartitionLeft( M& A, M& AL, M& AR, int widthAR )
{
#ifndef RELEASE
    PushCallStack("PartitionLeft [Matrix]");
    if( widthAR < 0 )
        throw std::logic_error("Width of right partition must be non-negative");
#endif
    widthAR = std::min(widthAR,A.Width());
    const int widthAL = A.Width()-widthAR;
    AL.View( A, 0, 0,       A.Height(), widthAL );
    AR.View( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
PartitionLeft( DM& A, DM& AL, DM& AR, int widthAR )
{
#ifndef RELEASE
    PushCallStack("PartitionLeft [DistMatrix]");
    if( widthAR < 0 )
        throw std::logic_error("Width of right partition must be non-negative");
#endif
    widthAR = std::min(widthAR,A.Width());
    const int widthAL = A.Width()-widthAR;
    AL.View( A, 0, 0,       A.Height(), widthAL );
    AR.View( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
PartitionRight( M& A, M& AL, M& AR, int widthAL )
{
#ifndef RELEASE
    PushCallStack("PartitionRight [Matrix]");
    if( widthAL < 0 )
        throw std::logic_error("Width of left partition must be non-negative");
#endif
    widthAL = std::min(widthAL,A.Width());
    const int widthAR = A.Width()-widthAL;
    AL.View( A, 0, 0,       A.Height(), widthAL );
    AR.View( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
PartitionRight( DM& A, DM& AL, DM& AR, int widthAL )
{
#ifndef RELEASE
    PushCallStack("PartitionRight [DistMatrix]");
    if( widthAL < 0 )
        throw std::logic_error("Width of left partition must be non-negative");
#endif
    widthAL = std::min(widthAL,A.Width());
    const int widthAR = A.Width()-widthAL;
    AL.View( A, 0, 0,       A.Height(), widthAL );
    AR.View( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
PartitionUpDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, int diagABR )
{
#ifndef RELEASE
    PushCallStack("PartitionUpDiagonal [Matrix]");
#endif
    PartitionUpLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, diagABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
PartitionUpDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, int diagABR )
{
#ifndef RELEASE
    PushCallStack("PartitionUpDiagonal [DistMatrix]");
#endif
    PartitionUpLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, diagABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
PartitionUpLeftDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, int diagABR )
{
#ifndef RELEASE
    PushCallStack("PartitionUpLeftDiagonal [Matrix]");
    if( diagABR < 0 )
        throw std::logic_error("Bottom-right size must be non-negative");
#endif
    const int minDim = std::min( A.Height(), A.Width() );
    diagABR = std::min(diagABR,minDim);
    const int sizeATL = minDim - diagABR;
    const int remHeight = A.Height()-sizeATL;
    const int remWidth = A.Width()-sizeATL;
    ATL.View( A, 0,       0,       sizeATL,   sizeATL  );
    ATR.View( A, 0,       sizeATL, sizeATL,   remWidth );
    ABL.View( A, sizeATL, 0,       remHeight, sizeATL  );
    ABR.View( A, sizeATL, sizeATL, remHeight, remWidth );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
PartitionUpLeftDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, int diagABR )
{
#ifndef RELEASE
    PushCallStack("PartitionUpLeftDiagonal [DistMatrix]");
    if( diagABR < 0 )
        throw std::logic_error("Bottom-right size must be non-negative");
#endif
    const int minDim = std::min( A.Height(), A.Width() );
    diagABR = std::min(diagABR,minDim);
    const int sizeATL = minDim - diagABR;
    const int remHeight = A.Height()-sizeATL;
    const int remWidth = A.Width()-sizeATL;
    ATL.View( A, 0,       0,       sizeATL,   sizeATL  );
    ATR.View( A, 0,       sizeATL, sizeATL,   remWidth );
    ABL.View( A, sizeATL, 0,       remHeight, sizeATL  );
    ABR.View( A, sizeATL, sizeATL, remHeight, remWidth );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
PartitionUpRightDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, int diagABR )
{
#ifndef RELEASE
    PushCallStack("PartitionUpRightDiagonal [Matrix]");
    if( diagABR < 0 )
        throw std::logic_error("Bottom-right size must be non-negative");
#endif
    const int minDim = std::min(A.Height(),A.Width());
    diagABR = std::min(diagABR,minDim);
    const int remHeight = A.Height()-diagABR;
    const int remWidth = A.Width()-diagABR;
    ATL.View( A, 0,         0,        remHeight, remWidth );
    ATR.View( A, 0,         remWidth, remHeight, diagABR  );
    ABL.View( A, remHeight, 0,        diagABR,   remWidth );
    ABR.View( A, remHeight, remWidth, diagABR,   diagABR  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
PartitionUpRightDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, int diagABR )
{
#ifndef RELEASE
    PushCallStack("PartitionUpRightDiagonal [DistMatrix]");
    if( diagABR < 0 )
        throw std::logic_error("Bottom-right size must be non-negative");
#endif
    const int minDim = std::min(A.Height(),A.Width());
    diagABR = std::min(diagABR,minDim);
    const int remHeight = A.Height()-diagABR;
    const int remWidth = A.Width()-diagABR;
    ATL.View( A, 0,         0,        remHeight, remWidth );
    ATR.View( A, 0,         remWidth, remHeight, diagABR  );
    ABL.View( A, remHeight, 0,        diagABR,   remWidth );
    ABR.View( A, remHeight, remWidth, diagABR,   diagABR  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
PartitionDownDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, int diagATL )
{
#ifndef RELEASE
    PushCallStack("PartitionDownDiagonal [Matrix]");
#endif
    PartitionDownLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, diagATL );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
PartitionDownDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, int diagATL )
{
#ifndef RELEASE
    PushCallStack("PartitionDownDiagonal [DistMatrix]");
#endif
    PartitionDownLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, diagATL );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
PartitionDownLeftDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, int diagATL )
{
#ifndef RELEASE
    PushCallStack("PartitionDownLeftDiagonal [Matrix]");
    if( diagATL < 0 )
        throw std::logic_error("Top-left size must be non-negative");
#endif
    const int minDim = std::min(A.Height(),A.Width());
    diagATL = std::min(diagATL,minDim);
    const int heightABR = A.Height()-diagATL;
    const int widthABR = A.Width()-diagATL;
    ATL.View( A, 0,       0,       diagATL,   diagATL  );
    ATR.View( A, 0,       diagATL, diagATL,   widthABR );
    ABL.View( A, diagATL, 0,       heightABR, diagATL  );
    ABR.View( A, diagATL, diagATL, heightABR, widthABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
PartitionDownLeftDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, int diagATL )
{
#ifndef RELEASE
    PushCallStack("PartitionDownLeftDiagonal [DistMatrix]");
    if( diagATL < 0 )
        throw std::logic_error("Top-left size must be non-negative");
#endif
    const int minDim = std::min(A.Height(),A.Width());
    diagATL = std::min(diagATL,minDim);
    const int heightABR = A.Height()-diagATL;
    const int widthABR = A.Width()-diagATL;
    ATL.View( A, 0,       0,       diagATL,   diagATL  );
    ATR.View( A, 0,       diagATL, diagATL,   widthABR );
    ABL.View( A, diagATL, 0,       heightABR, diagATL  );
    ABR.View( A, diagATL, diagATL, heightABR, widthABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
PartitionDownRightDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, int diagATL )
{
#ifndef RELEASE
    PushCallStack("PartitionDownRightDiagonal [Matrix]");
    if( diagATL < 0 )
        throw std::logic_error("Top-left size must be non-negative");
#endif
    const int minDim = std::min( A.Height(), A.Width() );
    diagATL = std::min(diagATL,minDim);
    const int sizeABR = minDim-diagATL;
    const int remHeight = A.Height()-sizeABR;
    const int remWidth = A.Width()-sizeABR;
    ATL.View( A, 0,         0,        remHeight, remWidth );
    ATR.View( A, 0,         remWidth, remHeight, sizeABR  );
    ABL.View( A, remHeight, 0,        sizeABR,   remWidth );
    ABR.View( A, remHeight, remWidth, sizeABR,   sizeABR  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
PartitionDownRightDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, int diagATL )
{
#ifndef RELEASE
    PushCallStack("PartitionDownRightDiagonal [DistMatrix]");
    if( diagATL < 0 )
        throw std::logic_error("Top-left size must be non-negative");
#endif
    const int minDim = std::min( A.Height(), A.Width() );
    diagATL = std::min(diagATL,minDim);
    const int sizeABR = minDim-diagATL;
    const int remHeight = A.Height()-sizeABR;
    const int remWidth = A.Width()-sizeABR;
    ATL.View( A, 0,         0,        remHeight, remWidth );
    ATR.View( A, 0,         remWidth, remHeight, sizeABR  );
    ABL.View( A, remHeight, 0,        sizeABR,   remWidth );
    ABR.View( A, remHeight, remWidth, sizeABR,   sizeABR  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LockedPartitionUp
( const M& A, M& AT,
              M& AB, int heightAB )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionUp [Matrix]");
    if( heightAB < 0 )
        throw std::logic_error
        ("Height of bottom partition must be non-negative");
#endif
    heightAB = std::min(heightAB,A.Height());
    const int heightAT = A.Height()-heightAB;
    AT.LockedView( A, 0,        0, heightAT, A.Width() );
    AB.LockedView( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
LockedPartitionUp
( const DM& A, DM& AT,
               DM& AB, int heightAB )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionUp [DistMatrix]");
    if( heightAB < 0 )
        throw std::logic_error
        ("Height of bottom partition must be non-negative");
#endif
    heightAB = std::min(heightAB,A.Height());
    const int heightAT = A.Height()-heightAB;
    AT.LockedView( A, 0,        0, heightAT, A.Width() );
    AB.LockedView( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LockedPartitionDown
( const M& A, M& AT,
              M& AB, int heightAT ) 
{
#ifndef RELEASE
    PushCallStack("LockedPartitionDown [Matrix]");
    if( heightAT < 0 )
        throw std::logic_error("Height of top partition must be non-negative");
#endif
    heightAT = std::min(heightAT,A.Height());
    const int heightAB = A.Height()-heightAT;
    AT.LockedView( A, 0,        0, heightAT, A.Width() );
    AB.LockedView( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
LockedPartitionDown
( const DM& A, DM& AT,
               DM& AB, int heightAT )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionDown [DistMatrix]");
    if( heightAT < 0 )
        throw std::logic_error("Height of top partition must be non-negative");
#endif
    heightAT = std::min(heightAT,A.Height());
    const int heightAB = A.Height()-heightAT;
    AT.LockedView( A, 0,        0, heightAT, A.Width() );
    AB.LockedView( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LockedPartitionLeft( const M& A, M& AL, M& AR, int widthAR )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionLeft [Matrix]");
    if( widthAR < 0 )
        throw std::logic_error("Width of right partition must be non-negative");
#endif
    widthAR = std::min(widthAR,A.Width());
    const int widthAL = A.Width()-widthAR;
    AL.LockedView( A, 0, 0,       A.Height(), widthAL );
    AR.LockedView( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
LockedPartitionLeft( const DM& A, DM& AL, DM& AR, int widthAR )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionLeft [DistMatrix]");
    if( widthAR < 0 )
        throw std::logic_error("Width of right partition must be non-negative");
#endif
    widthAR = std::min(widthAR,A.Width());
    const int widthAL = A.Width()-widthAR;
    AL.LockedView( A, 0, 0,       A.Height(), widthAL );
    AR.LockedView( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LockedPartitionRight( const M& A, M& AL, M& AR, int widthAL )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionRight [Matrix]");
    if( widthAL < 0 )
        throw std::logic_error("Width of left partition must be non-negative");
#endif
    widthAL = std::min(widthAL,A.Width());
    const int widthAR = A.Width()-widthAL;
    AL.LockedView( A, 0, 0,       A.Height(), widthAL );
    AR.LockedView( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
LockedPartitionRight( const DM& A, DM& AL, DM& AR, int widthAL )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionRight [DistMatrix]");
    if( widthAL < 0 )
        throw std::logic_error("Width of left partition must be non-negative");
#endif
    widthAL = std::min(widthAL,A.Width());
    const int widthAR = A.Width()-widthAL;
    AL.LockedView( A, 0, 0,       A.Height(), widthAL );
    AR.LockedView( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LockedPartitionUpDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, int diagABR )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionUpDiagonal [Matrix]");
#endif
    LockedPartitionUpLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, diagABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
LockedPartitionUpDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, int diagABR )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionUpDiagonal [DistMatrix]");
#endif
    LockedPartitionUpLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, diagABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LockedPartitionUpLeftDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, int diagABR )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionUpLeftDiagonal [Matrix]");
    if( diagABR < 0 )
        throw std::logic_error("Bottom-right size must be non-negative");
#endif
    const int minDim = std::min( A.Height(), A.Width() );
    diagABR = std::min(diagABR,minDim);
    const int sizeATL = minDim - diagABR;
    const int remHeight = A.Height()-sizeATL;
    const int remWidth = A.Width()-sizeATL;
    ATL.LockedView( A, 0,       0,       sizeATL,   sizeATL  );
    ATR.LockedView( A, 0,       sizeATL, sizeATL,   remWidth );
    ABL.LockedView( A, sizeATL, 0,       remHeight, sizeATL  );
    ABR.LockedView( A, sizeATL, sizeATL, remHeight, remWidth );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
LockedPartitionUpLeftDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, int diagABR )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionUpLeftDiagonal [DistMatrix]");
    if( diagABR < 0 )
        throw std::logic_error("Bottom-right size must be non-negative");
#endif
    const int minDim = std::min( A.Height(), A.Width() );
    diagABR = std::min(diagABR,minDim);
    const int sizeATL = minDim - diagABR;
    const int remHeight = A.Height()-sizeATL;
    const int remWidth = A.Width()-sizeATL;
    ATL.LockedView( A, 0,       0,       sizeATL,   sizeATL  );
    ATR.LockedView( A, 0,       sizeATL, sizeATL,   remWidth );
    ABL.LockedView( A, sizeATL, 0,       remHeight, sizeATL  );
    ABR.LockedView( A, sizeATL, sizeATL, remHeight, remWidth );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LockedPartitionUpRightDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, int diagABR )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionUpRightDiagonal [Matrix]");
    if( diagABR < 0 )
        throw std::logic_error("Bottom-right size must be non-negative");
#endif
    const int minDim = std::min(A.Height(),A.Width());
    diagABR = std::min(diagABR,minDim);
    const int remHeight = A.Height()-diagABR;
    const int remWidth = A.Width()-diagABR;
    ATL.View( A, 0,         0,        remHeight, remWidth );
    ATR.View( A, 0,         remWidth, remHeight, diagABR  );
    ABL.View( A, remHeight, 0,        diagABR,   remWidth );
    ABR.View( A, remHeight, remWidth, diagABR,   diagABR  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
LockedPartitionUpRightDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, int diagABR )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionUpRightDiagonal [DistMatrix]");
    if( diagABR < 0 )
        throw std::logic_error("Bottom-right size must be non-negative");
#endif
    const int minDim = std::min(A.Height(),A.Width());
    diagABR = std::min(diagABR,minDim);
    const int remHeight = A.Height()-diagABR;
    const int remWidth = A.Width()-diagABR;
    ATL.View( A, 0,         0,        remHeight, remWidth );
    ATR.View( A, 0,         remWidth, remHeight, diagABR  );
    ABL.View( A, remHeight, 0,        diagABR,   remWidth );
    ABR.View( A, remHeight, remWidth, diagABR,   diagABR  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LockedPartitionDownDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, int diagATL )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionDownDiagonal [Matrix]");
#endif
    LockedPartitionDownLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, diagATL );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
LockedPartitionDownDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, int diagATL )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionDownDiagonal [DistMatrix]");
#endif
    LockedPartitionDownLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, diagATL );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LockedPartitionDownLeftDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, int diagATL )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionDownLeftDiagonal [Matrix]");
    if( diagATL < 0 )
        throw std::logic_error("Top-left size must be non-negative");
#endif
    const int minDim = std::min(A.Height(),A.Width());
    diagATL = std::min(diagATL,minDim);
    const int heightABR = A.Height()-diagATL;
    const int widthABR = A.Width()-diagATL;
    ATL.LockedView( A, 0,       0,       diagATL,   diagATL  );
    ATR.LockedView( A, 0,       diagATL, diagATL,   widthABR );
    ABL.LockedView( A, diagATL, 0,       heightABR, diagATL  );
    ABR.LockedView( A, diagATL, diagATL, heightABR, widthABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
LockedPartitionDownLeftDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, int diagATL )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionDownLeftDiagonal [DistMatrix]");
    if( diagATL < 0 )
        throw std::logic_error("Top-left size must be non-negative");
#endif
    const int minDim = std::min(A.Height(),A.Width());
    diagATL = std::min(diagATL,minDim);
    const int heightABR = A.Height()-diagATL;
    const int widthABR = A.Width()-diagATL;
    ATL.LockedView( A, 0,       0,       diagATL,   diagATL  );
    ATR.LockedView( A, 0,       diagATL, diagATL,   widthABR );
    ABL.LockedView( A, diagATL, 0,       heightABR, diagATL  );
    ABR.LockedView( A, diagATL, diagATL, heightABR, widthABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LockedPartitionDownRightDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, int diagATL )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionDownRightDiagonal [Matrix]");
    if( diagATL < 0 )
        throw std::logic_error("Top-left size must be non-negative");
#endif
    const int minDim = std::min( A.Height(), A.Width() );
    diagATL = std::min(diagATL,minDim);
    const int sizeABR = minDim-diagATL;
    const int remHeight = A.Height()-sizeABR;
    const int remWidth = A.Width()-sizeABR;
    ATL.View( A, 0,         0,        remHeight, remWidth );
    ATR.View( A, 0,         remWidth, remHeight, sizeABR  );
    ABL.View( A, remHeight, 0,        sizeABR,   remWidth );
    ABR.View( A, remHeight, remWidth, sizeABR,   sizeABR  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
LockedPartitionDownRightDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, int diagATL )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionDownRightDiagonal [DistMatrix]");
    if( diagATL < 0 )
        throw std::logic_error("Top-left size must be non-negative");
#endif
    const int minDim = std::min( A.Height(), A.Width() );
    diagATL = std::min(diagATL,minDim);
    const int sizeABR = minDim-diagATL;
    const int remHeight = A.Height()-sizeABR;
    const int remWidth = A.Width()-sizeABR;
    ATL.View( A, 0,         0,        remHeight, remWidth );
    ATR.View( A, 0,         remWidth, remHeight, sizeABR  );
    ABL.View( A, remHeight, 0,        sizeABR,   remWidth );
    ABR.View( A, remHeight, remWidth, sizeABR,   sizeABR  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
RepartitionUp
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2, int A1Height )
{
#ifndef RELEASE
    PushCallStack("RepartitionUp [Matrix]");
    if( (AT.Buffer() + AT.Height()) != AB.Buffer() )
        throw std::logic_error("Noncontiguous 2x1 array of matrices");
#endif
    A1Height = std::min(A1Height,AT.Height());
    const int offset = AT.Height()-A1Height; 
    A0.View( AT, 0,      0, offset,   AT.Width() );
    A1.View( AT, offset, 0, A1Height, AT.Width() );
    A2.View( AB );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
RepartitionUp
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2, int A1Height )
{
#ifndef RELEASE
    PushCallStack("RepartitionUp [DistMatrix]");
    if( (AT.LocalMatrix().Buffer() + AT.LocalHeight()) != 
         AB.LocalMatrix().Buffer() )
        throw std::logic_error
        ("Noncontiguous 2x1 array of distributed matrices");
#endif
    A1Height = std::min(A1Height,AT.Height());
    const int offset = AT.Height()-A1Height; 
    A0.View( AT, 0,      0, offset,   AT.Width() );
    A1.View( AT, offset, 0, A1Height, AT.Width() );
    A2.View( AB );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
RepartitionDown
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2, int A1Height )
{
#ifndef RELEASE
    PushCallStack("RepartitionDown [Matrix]");
    if( (AT.Buffer() + AT.Height()) != AB.Buffer() )
        throw std::logic_error("Noncontiguous 2x1 array of matrices");
#endif
    A1Height = std::min(A1Height,AB.Height());
    const int offset = AB.Height()-A1Height; 
    A0.View( AT );
    A1.View( AB, 0,        0, A1Height, AB.Width() );
    A2.View( AB, A1Height, 0, offset,   AB.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
RepartitionDown
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2, int A1Height )
{
#ifndef RELEASE
    PushCallStack("RepartitionDown [DistMatrix]");
    if( (AT.LocalMatrix().Buffer() + AT.LocalHeight()) != 
         AB.LocalMatrix().Buffer() )
        throw std::logic_error
        ("Noncontiguous 2x1 array of distributed matrices");
#endif
    A1Height = std::min(A1Height,AB.Height());
    const int offset = AB.Height()-A1Height; 
    A0.View( AT );
    A1.View( AB, 0,        0, A1Height, AB.Width() );
    A2.View( AB, A1Height, 0, offset, AB.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
RepartitionLeft
( M& AL, M& AR,
  M& A0, M& A1, M& A2, int A1Width )
{
#ifndef RELEASE
    PushCallStack("RepartitionLeft [Matrix]");
    if( (AL.Buffer() + AL.Width()*AL.LDim()) != AR.Buffer() )
        throw std::logic_error("Noncontiguous 1x2 array of matrices");
#endif
    A1Width = std::min(A1Width,AL.Width());
    const int offset = AL.Width()-A1Width;
    A0.View( AL, 0, 0,      AL.Height(), offset   );
    A1.View( AL, 0, offset, AL.Height(), A1Width  );
    A2.View( AR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
RepartitionLeft
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2, int A1Width )
{
#ifndef RELEASE
    PushCallStack("RepartitionLeft [DistMatrix]");
    if( (AL.LocalMatrix().Buffer() + AL.LocalWidth()*AL.LocalLDim())
         != AR.LocalMatrix().Buffer() )
        throw std::logic_error
        ("Noncontiguous 1x2 array of distributed matrices");
#endif
    A1Width = std::min(A1Width,AL.Width());
    const int offset = AL.Width()-A1Width;
    A0.View( AL, 0, 0,      AL.Height(), offset  );
    A1.View( AL, 0, offset, AL.Height(), A1Width );
    A2.View( AR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
RepartitionRight
( M& AL, M& AR,
  M& A0, M& A1, M& A2, int A1Width )
{
#ifndef RELEASE
    PushCallStack("RepartitionRight [Matrix]");
    if( (AL.Buffer() + AL.Width()*AL.LDim()) != AR.Buffer() )
        throw std::logic_error("Noncontiguous 1x2 array of matrices");
#endif
    A1Width = std::min(A1Width,AR.Width());
    const int offset = AR.Width()-A1Width;
    A0.View( AL );
    A1.View( AR, 0, 0,       AR.Height(), A1Width );
    A2.View( AR, 0, A1Width, AR.Height(), offset  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
RepartitionRight
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2, int A1Width )
{
#ifndef RELEASE
    PushCallStack("RepartitionRight [DistMatrix]");
    if( (AL.LocalMatrix().Buffer() + AL.LocalWidth()*AL.LocalLDim()) 
         != AR.LocalMatrix().Buffer() )
        throw std::logic_error
        ("Noncontiguous 1x2 array of distributed matrices");
#endif
    A1Width = std::min(A1Width,AR.Width());
    const int offset = AR.Width()-A1Width;
    A0.View( AL );
    A1.View( AR, 0, 0,       AR.Height(), A1Width );
    A2.View( AR, 0, A1Width, AR.Height(), offset  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
RepartitionUpDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22, int bsize )
{
#ifndef RELEASE
    PushCallStack("RepartitionUpDiagonal [Matrix]");
    if( (ATL.Buffer() + ATL.Height()) != ABL.Buffer() ||
        (ATR.Buffer() + ATR.Height()) != ABR.Buffer() ||
        (ATL.Buffer() + ATL.Width()*ATL.LDim()) != ATR.Buffer() ||
        (ABL.Buffer() + ABL.Width()*ABL.LDim()) != ABR.Buffer()    )
        throw std::logic_error("Noncontiguous 2x2 grid of matrices");
#endif
    bsize = std::min(bsize,std::min(ATL.Height(),ATL.Width()));
    const int vOffset = ATL.Height()-bsize;
    const int hOffset = ATL.Width()-bsize;
    A00.View( ATL, 0,       0,       vOffset,      hOffset     );
    A01.View( ATL, 0,       hOffset, vOffset,      bsize       );
    A02.View( ATR, 0,       0,       vOffset,      ATR.Width() );
    A10.View( ATL, vOffset, 0,       bsize,        hOffset     );
    A11.View( ATL, vOffset, hOffset, bsize,        bsize       );
    A12.View( ATR, vOffset, 0,       bsize,        ATR.Width() );
    A20.View( ABL, 0,       0,       ABL.Height(), hOffset     );
    A21.View( ABL, 0,       hOffset, ABL.Height(), bsize       );
    A22.View( ABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
RepartitionUpDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22, int bsize )
{
#ifndef RELEASE
    PushCallStack("RepartitionUpDiagonal [DistMatrix]");
    if( (ATL.LocalMatrix().Buffer() + ATL.LocalHeight()) != 
         ABL.LocalMatrix().Buffer() ||
        (ATR.LocalMatrix().Buffer() + ATR.LocalHeight()) != 
         ABR.LocalMatrix().Buffer() ||
        (ATL.LocalMatrix().Buffer() + ATL.LocalWidth()*ATL.LocalLDim())
         != ATR.LocalMatrix().Buffer() ||
        (ABL.LocalMatrix().Buffer() + ABL.LocalWidth()*ABL.LocalLDim())
         != ABR.LocalMatrix().Buffer() )
        throw std::logic_error
        ("Noncontiguous 2x2 grid of distributed matrices");
#endif
    bsize = std::min(bsize,std::min(ATL.Height(),ATL.Width()));
    const int vOffset = ATL.Height()-bsize;
    const int hOffset = ATL.Width()-bsize;
    A00.View( ATL, 0,       0,       vOffset,      hOffset     );
    A01.View( ATL, 0,       hOffset, vOffset,      bsize       );
    A02.View( ATR, 0,       0,       vOffset,      ATR.Width() );
    A10.View( ATL, vOffset, 0,       bsize,        hOffset     );
    A11.View( ATL, vOffset, hOffset, bsize,        bsize       );
    A12.View( ATR, vOffset, 0,       bsize,        ATR.Width() );
    A20.View( ABL, 0,       0,       ABL.Height(), hOffset     );
    A21.View( ABL, 0,       hOffset, ABL.Height(), bsize       );
    A22.View( ABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
RepartitionDownDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22, int bsize )
{
#ifndef RELEASE
    PushCallStack("RepartitionDownDiagonal [Matrix]");
    if( (ATL.Buffer() + ATL.Height()) != ABL.Buffer() ||
        (ATR.Buffer() + ATR.Height()) != ABR.Buffer() ||
        (ATL.Buffer() + ATL.Width()*ATL.LDim()) != ATR.Buffer() ||
        (ABL.Buffer() + ABL.Width()*ABL.LDim()) != ABR.Buffer()    )
        throw std::logic_error("Noncontiguous 2x2 grid of matrices");
#endif
    bsize = std::min(bsize,std::min(ABR.Height(),ABR.Width()));
    const int vOffset = ABR.Height()-bsize;
    const int hOffset = ABR.Width()-bsize;
    A00.View( ATL );
    A01.View( ATR, 0,     0,     ATL.Height(), bsize       );
    A02.View( ATR, 0,     bsize, ATL.Height(), hOffset     );
    A10.View( ABL, 0,     0,     bsize,        ABL.Width() );
    A11.View( ABR, 0,     0,     bsize,        bsize       );
    A12.View( ABR, 0,     bsize, bsize,        hOffset     );
    A20.View( ABL, bsize, 0,     vOffset,      ABL.Width() );
    A21.View( ABR, bsize, 0,     vOffset,      bsize       );
    A22.View( ABR, bsize, bsize, vOffset,      hOffset     );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
RepartitionDownDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22, int bsize )
{
#ifndef RELEASE
    PushCallStack("RepartitionDownDiagonal [DistMatrix]");
    if( (ATL.LocalMatrix().Buffer() + ATL.LocalHeight()) != 
         ABL.LocalMatrix().Buffer() ||
        (ATR.LocalMatrix().Buffer() + ATR.LocalHeight()) != 
         ABR.LocalMatrix().Buffer() ||
        (ATL.LocalMatrix().Buffer() + ATL.LocalWidth()*ATL.LocalLDim()) !=
         ATR.LocalMatrix().Buffer() ||
        (ABL.LocalMatrix().Buffer() + ABL.LocalWidth()*ABL.LocalLDim()) != 
         ABR.LocalMatrix().Buffer() )
        throw std::logic_error
        ("Noncontiguous 2x2 grid of distributed matrices");
#endif
    bsize = std::min(bsize,std::min(ABR.Height(),ABR.Width()));
    const int vOffset = ABR.Height()-bsize;
    const int hOffset = ABR.Width()-bsize;
    A00.View( ATL );
    A01.View( ATR, 0,     0,     ATL.Height(), bsize       );
    A02.View( ATR, 0,     bsize, ATL.Height(), hOffset     );
    A10.View( ABL, 0,     0,     bsize,        ABL.Width() );
    A11.View( ABR, 0,     0,     bsize,        bsize       );
    A12.View( ABR, 0,     bsize, bsize,        hOffset     );
    A20.View( ABL, bsize, 0,     vOffset,      ABL.Width() );
    A21.View( ABR, bsize, 0,     vOffset,      bsize       );
    A22.View( ABR, bsize, bsize, vOffset,      hOffset     );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LockedRepartitionUp
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2, int A1Height )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionUp [Matrix]");
    if( (AT.LockedBuffer() + AT.Height()) != AB.LockedBuffer() )
        throw std::logic_error("Noncontiguous 2x1 array of matrices");
#endif
    A1Height = std::min(A1Height,AT.Height());
    const int offset = AT.Height()-A1Height;
    A0.LockedView( AT, 0,      0, offset,   AT.Width() );
    A1.LockedView( AT, offset, 0, A1Height, AT.Width() );
    A2.LockedView( AB );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
LockedRepartitionUp
( const DM& AT, DM& A0,
                DM& A1,
  const DM& AB, DM& A2, int A1Height )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionUp [DistMatrix]");
    if( (AT.LockedLocalMatrix().LockedBuffer() + AT.LocalHeight()) != 
         AB.LockedLocalMatrix().LockedBuffer() )
        throw std::logic_error
        ("Noncontiguous 2x1 array of distributed matrices");
#endif
    A1Height = std::min(A1Height,AT.Height());
    const int offset = AT.Height()-A1Height;
    A0.LockedView( AT, 0,      0, offset,   AT.Width() );
    A1.LockedView( AT, offset, 0, A1Height, AT.Width() );
    A2.LockedView( AB );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LockedRepartitionDown
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2, int A1Height )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionDown [Matrix]");
    if( (AT.LockedBuffer() + AT.Height()) != AB.LockedBuffer() )
        throw std::logic_error("Noncontiguous 2x1 array of matrices");
#endif
    A1Height = std::min(A1Height,AB.Height());
    const int offset = AB.Height()-A1Height;
    A0.LockedView( AT );
    A1.LockedView( AB, 0,        0, A1Height, AB.Width() );
    A2.LockedView( AB, A1Height, 0, offset,   AB.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
LockedRepartitionDown
( const DM& AT, DM& A0,
                DM& A1,
  const DM& AB, DM& A2, int A1Height )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionDown [DistMatrix]");
    if( (AT.LockedLocalMatrix().LockedBuffer() + AT.LocalHeight()) != 
         AB.LockedLocalMatrix().LockedBuffer() )
        throw std::logic_error
        ("Noncontiguous 2x1 array of distributed matrices");
#endif
    A1Height = std::min(A1Height,AB.Height());
    const int offset = AB.Height()-A1Height;
    A0.LockedView( AT );
    A1.LockedView( AB, 0,        0, A1Height, AB.Width() );
    A2.LockedView( AB, A1Height, 0, offset,   AB.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LockedRepartitionLeft
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2, int A1Width )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionLeft [Matrix]");
    if( (AL.LockedBuffer() + AL.Width()*AL.LDim()) != AR.LockedBuffer() )
        throw std::logic_error("Noncontiguous 1x2 array of matrices");
#endif
    A1Width = std::min(A1Width,AL.Width());
    const int offset = AL.Width()-A1Width;
    A0.LockedView( AL, 0, 0,      AL.Height(), offset  );
    A1.LockedView( AL, 0, offset, AL.Height(), A1Width );
    A2.LockedView( AR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
LockedRepartitionLeft
( const DM& AL, const DM& AR,
  DM& A0, DM& A1, DM& A2, int A1Width )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionLeft [DistMatrix]");
    if( (AL.LockedLocalMatrix().LockedBuffer() + AL.LocalWidth()*AL.LocalLDim())
         != AR.LockedLocalMatrix().LockedBuffer() )
        throw std::logic_error
        ("Noncontiguous 1x1 array of distributed matrices");
#endif
    A1Width = std::min(A1Width,AL.Width());
    const int offset = AL.Width()-A1Width;
    A0.LockedView( AL, 0, 0,      AL.Height(), offset  );
    A1.LockedView( AL, 0, offset, AL.Height(), A1Width );
    A2.LockedView( AR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LockedRepartitionRight
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2, int A1Width )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionRight [Matrix]");
    if( (AL.LockedBuffer() + AL.Width()*AL.LDim()) != AR.LockedBuffer() )
        throw std::logic_error("Noncontiguous 1x2 array of matrices");
#endif
    A1Width = std::min(A1Width,AR.Width());
    const int offset = AR.Width()-A1Width;
    A0.LockedView( AL );
    A1.LockedView( AR, 0, 0,       AR.Height(), A1Width );
    A2.LockedView( AR, 0, A1Width, AR.Height(), offset  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
LockedRepartitionRight
( const DM& AL, const DM& AR,
  DM& A0, DM& A1, DM& A2, int A1Width )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionRight [DistMatrix]");
    if( (AL.LockedLocalMatrix().LockedBuffer() + AL.LocalWidth()*AL.LocalLDim())
         != AR.LockedLocalMatrix().LockedBuffer() )
        throw std::logic_error
        ("Noncontiguous 1x2 DistMatrices in LockedRepartitionRight");
#endif
    A1Width = std::min(A1Width,AR.Width());
    const int offset = AR.Width()-A1Width;
    A0.LockedView( AL );
    A1.LockedView( AR, 0, 0,       AR.Height(), A1Width );
    A2.LockedView( AR, 0, A1Width, AR.Height(), offset  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LockedRepartitionUpDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22, int bsize )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionUpDiagonal [Matrix]");
    if( (ATL.LockedBuffer() + ATL.Height()) != ABL.LockedBuffer() ||
        (ATR.LockedBuffer() + ATR.Height()) != ABR.LockedBuffer() ||
        (ATL.LockedBuffer() + ATL.Width()*ATL.LDim()) != ATR.LockedBuffer() ||
        (ABL.LockedBuffer() + ABL.Width()*ABL.LDim()) != ABR.LockedBuffer() )
        throw std::logic_error("Noncontiguous 2x2 grid of matrices");
#endif
    bsize = std::min(bsize,std::min(ATL.Height(),ATL.Width()));
    const int vOffset = ATL.Height()-bsize;
    const int hOffset = ATL.Width()-bsize;
    A00.LockedView( ATL, 0,       0,       vOffset,      hOffset     );
    A01.LockedView( ATL, 0,       hOffset, vOffset,      bsize       );
    A02.LockedView( ATR, 0,       0,       vOffset,      ATR.Width() );
    A10.LockedView( ATL, vOffset, 0,       bsize,        hOffset     );
    A11.LockedView( ATL, vOffset, hOffset, bsize,        bsize       );
    A12.LockedView( ATR, vOffset, 0,       bsize,        ATR.Width() );
    A20.LockedView( ABL, 0,       0,       ABL.Height(), hOffset     );
    A21.LockedView( ABL, 0,       hOffset, ABL.Height(), bsize       );
    A22.LockedView( ABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
LockedRepartitionUpDiagonal
( const DM& ATL, const DM& ATR, DM& A00, DM& A01, DM& A02,
                                DM& A10, DM& A11, DM& A12,
  const DM& ABL, const DM& ABR, DM& A20, DM& A21, DM& A22, int bsize )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionUpDiagonal [DistMatrix]");
    if( (ATL.LockedLocalMatrix().LockedBuffer()+ATL.LocalHeight()) != 
         ABL.LockedLocalMatrix().LockedBuffer() ||
        (ATR.LockedLocalMatrix().LockedBuffer()+ATR.LocalHeight()) != 
         ABR.LockedLocalMatrix().LockedBuffer() ||
        (ATL.LockedLocalMatrix().LockedBuffer()+
         ATL.LocalWidth()*ATL.LocalLDim()) != 
         ATR.LockedLocalMatrix().LockedBuffer() ||
        (ABL.LockedLocalMatrix().LockedBuffer()+
         ABL.LocalWidth()*ABL.LocalLDim()) !=
         ABR.LockedLocalMatrix().LockedBuffer() )
        throw std::logic_error
        ("Noncontiguous 2x2 grid of distributed matrices");
#endif
    bsize = std::min(bsize,std::min(ATL.Height(),ATL.Width()));
    const int vOffset = ATL.Height()-bsize;
    const int hOffset = ATL.Width()-bsize;
    A00.LockedView( ATL, 0,       0,       vOffset,      hOffset     );
    A01.LockedView( ATL, 0,       hOffset, vOffset,      bsize       );
    A02.LockedView( ATR, 0,       0,       vOffset,      ATR.Width() );
    A10.LockedView( ATL, vOffset, 0,       bsize,        hOffset     );
    A11.LockedView( ATL, vOffset, hOffset, bsize,        bsize       );
    A12.LockedView( ATR, vOffset, 0,       bsize,        ATR.Width() );
    A20.LockedView( ABL, 0,       0,       ABL.Height(), hOffset     );
    A21.LockedView( ABL, 0,       hOffset, ABL.Height(), bsize       );
    A22.LockedView( ABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LockedRepartitionDownDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22, int bsize )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionDownDiagonal [Matrix]");
    if( (ATL.LockedBuffer() + ATL.Height()) != ABL.LockedBuffer() ||
        (ATR.LockedBuffer() + ATR.Height()) != ABR.LockedBuffer() ||
        (ATL.LockedBuffer() + ATL.Width()*ATL.LDim()) != ATR.LockedBuffer() ||
        (ABL.LockedBuffer() + ABL.Width()*ABL.LDim()) != ABR.LockedBuffer() )
        throw std::logic_error("Noncontiguous 2x2 grid of matrices");
#endif
    bsize = std::min(bsize,std::min(ABR.Height(),ABR.Width()));
    const int vOffset = ABR.Height()-bsize;
    const int hOffset = ABR.Width()-bsize;
    A00.LockedView( ATL );
    A01.LockedView( ATR, 0,     0,     ATL.Height(), bsize       );
    A02.LockedView( ATR, 0,     bsize, ATL.Height(), hOffset     ); 
    A10.LockedView( ABL, 0,     0,     bsize,        ABL.Width() );
    A11.LockedView( ABR, 0,     0,     bsize,        bsize       );
    A12.LockedView( ABR, 0,     bsize, bsize,        hOffset     );
    A20.LockedView( ABL, bsize, 0,     vOffset,      ABL.Width() );
    A21.LockedView( ABR, bsize, 0,     vOffset,      bsize       );
    A22.LockedView( ABR, bsize, bsize, vOffset,      hOffset     );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
LockedRepartitionDownDiagonal
( const DM& ATL, const DM& ATR, DM& A00, DM& A01, DM& A02,
                                DM& A10, DM& A11, DM& A12,
  const DM& ABL, const DM& ABR, DM& A20, DM& A21, DM& A22, int bsize )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionDownDiagonal [DistMatrix]");
    if( (ATL.LockedLocalMatrix().LockedBuffer()+ATL.LocalHeight()) != 
         ABL.LockedLocalMatrix().LockedBuffer() ||
        (ATR.LockedLocalMatrix().LockedBuffer()+ATR.LocalHeight()) != 
         ABR.LockedLocalMatrix().LockedBuffer() ||
        (ATL.LockedLocalMatrix().LockedBuffer()+
         ATL.LocalWidth()*ATL.LocalLDim()) !=
         ATR.LockedLocalMatrix().LockedBuffer() ||
        (ABL.LockedLocalMatrix().LockedBuffer()+
         ABL.LocalWidth()*ABL.LocalLDim()) !=
         ABR.LockedLocalMatrix().LockedBuffer() )
        throw std::logic_error
        ("Noncontiguous 2x2 grid of distributed matrices");
#endif
    bsize = std::min(bsize,std::min(ABR.Height(),ABR.Width()));
    const int vOffset = ABR.Height()-bsize;
    const int hOffset = ABR.Width()-bsize;
    A00.LockedView( ATL );
    A01.LockedView( ATR, 0,     0,     ATR.Height(), bsize  );
    A02.LockedView( ATR, 0,     bsize, ATR.Height(), hOffset     );
    A10.LockedView( ABL, 0,     0,     bsize,        ABL.Width() );
    A11.LockedView( ABR, 0,     0,     bsize,        bsize       );
    A12.LockedView( ABR, 0,     bsize, bsize,        hOffset     );
    A20.LockedView( ABL, bsize, 0,     vOffset,      ABL.Width() );
    A21.LockedView( ABR, bsize, 0,     vOffset,      bsize       );
    A22.LockedView( ABR, bsize, bsize, vOffset,      hOffset     );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
SlidePartitionUp
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionUp [Matrix]");
#endif
    AT.View( A0 );
    AB.View2x1( A1,
                A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
SlidePartitionUp
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionUp [DistMatrix]");
#endif
    AT.View( A0 );
    AB.View2x1( A1,
                A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
SlidePartitionDown
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionDown [Matrix]");
#endif
    AT.View2x1( A0,
                A1 );
    AB.View( A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
SlidePartitionDown
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionDown [DistMatrix]");
#endif
    AT.View2x1( A0,
                A1 );
    AB.View( A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
SlidePartitionLeft
( M& AL, M& AR,
  M& A0, M& A1, M& A2 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionLeft [Matrix]");
#endif
    AL.View( A0 );
    AR.View1x2( A1, A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
SlidePartitionLeft
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionLeft [DistMatrix]");
#endif
    AL.View( A0 );
    AR.View1x2( A1, A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
SlidePartitionRight
( M& AL, M& AR,
  M& A0, M& A1, M& A2 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionRight [Matrix]");
#endif
    AL.View1x2( A0, A1 );
    AR.View( A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
SlidePartitionRight
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionRight [DistMatrix]");
#endif
    AL.View1x2( A0, A1 );
    AR.View( A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
SlidePartitionUpDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionUpDiagonal [Matrix]");
#endif
    ATL.View( A00 );
    ATR.View1x2( A01, A02 );
    ABL.View2x1( A10,
                 A20 );
    ABR.View2x2( A11, A12,
                 A21, A22 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
SlidePartitionUpDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionUpDiagonal [DistMatrix]");
#endif
    ATL.View( A00 );
    ATR.View1x2( A01, A02 );
    ABL.View2x1( A10,
                 A20 );
    ABR.View2x2( A11, A12,
                 A21, A22 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
SlidePartitionDownDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionDownDiagonal [Matrix]");
#endif
    ATL.View2x2( A00, A01,
                 A10, A11 );
    ATR.View2x1( A02,
                 A12 );
    ABL.View1x2( A20, A21 );
    ABR.View( A22 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
SlidePartitionDownDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionDownDiagonal [DistMatrix]");
#endif
    ATL.View2x2( A00, A01,
                 A10, A11 );
    ATR.View2x1( A02,
                 A12 );
    ABL.View1x2( A20, A21 );
    ABR.View( A22 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
SlideLockedPartitionUp
( M& AT, const M& A0,
         const M& A1,
  M& AB, const M& A2 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionUp [Matrix]");
#endif
    AT.LockedView( A0 );
    AB.LockedView2x1( A1,
                      A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
SlideLockedPartitionUp
( DM& AT, const DM& A0,
          const DM& A1,
  DM& AB, const DM& A2 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionUp [DistMatrix]");
#endif
    AT.LockedView( A0 );
    AB.LockedView2x1( A1,
                      A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
SlideLockedPartitionDown
( M& AT, const M& A0,
         const M& A1,
  M& AB, const M& A2 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionDown [Matrix]");
#endif
    AT.LockedView2x1( A0,
                      A1 );
    AB.LockedView( A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
SlideLockedPartitionDown
( DM& AT, const DM& A0,
          const DM& A1,
  DM& AB, const DM& A2 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionDown [DistMatrix]");
#endif
    AT.LockedView2x1( A0,
                      A1 );
    AB.LockedView( A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
SlideLockedPartitionLeft
( M& AL, M& AR,
  const M& A0, const M& A1, const M& A2 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionLeft [Matrix]");
#endif
    AL.LockedView( A0 );
    AR.LockedView1x2( A1, A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
SlideLockedPartitionLeft
( DM& AL, DM& AR,
  const DM& A0, const DM& A1, const DM& A2 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionLeft [DistMatrix]");
#endif
    AL.LockedView( A0 );
    AR.LockedView1x2( A1, A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
SlideLockedPartitionRight
( M& AL, M& AR,
  const M& A0, const M& A1, const M& A2 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionRight [Matrix]");
#endif
    AL.LockedView1x2( A0, A1 );
    AR.LockedView( A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
SlideLockedPartitionRight
( DM& AL, DM& AR,
  const DM& A0, const DM& A1, const DM& A2 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionRight [DistMatrix]");
#endif
    AL.LockedView1x2( A0, A1 );
    AR.LockedView( A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
SlideLockedPartitionUpDiagonal
( M& ATL, M& ATR, const M& A00, const M& A01, const M& A02,
                  const M& A10, const M& A11, const M& A12,
  M& ABL, M& ABR, const M& A20, const M& A21, const M& A22 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionUpDiagonal [Matrix]");
#endif
    ATL.LockedView( A00 );
    ATR.LockedView1x2( A01, A02 );
    ABL.LockedView2x1( A10,
                       A20 );
    ABR.LockedView2x2( A11, A12,
                       A21, A22 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
SlideLockedPartitionUpDiagonal
( DM& ATL, DM& ATR, const DM& A00, const DM& A01, const DM& A02,
                    const DM& A10, const DM& A11, const DM& A12,
  DM& ABL, DM& ABR, const DM& A20, const DM& A21, const DM& A22 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionUpDiagonal [DistMatrix]");
#endif
    ATL.LockedView( A00 );
    ATR.LockedView1x2( A01, A02 );
    ABL.LockedView2x1( A10,
                       A20 );
    ABR.LockedView2x2( A11, A12,
                       A21, A22 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
SlideLockedPartitionDownDiagonal
( M& ATL, M& ATR, const M& A00, const M& A01, const M& A02,
                  const M& A10, const M& A11, const M& A12,
  M& ABL, M& ABR, const M& A20, const M& A21, const M& A22 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionDownDiagonal [Matrix]");
#endif
    ATL.LockedView2x2( A00, A01,
                       A10, A11 );
    ATR.LockedView2x1( A02,
                       A12 );
    ABL.LockedView1x2( A20, A21 );
    ABR.LockedView( A22 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V>
inline void
SlideLockedPartitionDownDiagonal
( DM& ATL, DM& ATR, const DM& A00, const DM& A01, const DM& A02,
                    const DM& A10, const DM& A11, const DM& A12,
  DM& ABL, DM& ABR, const DM& A20, const DM& A21, const DM& A22 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionDownDiagonal [DistMatrix]");
#endif
    ATL.LockedView2x2( A00, A01,
                       A10, A11 );
    ATR.LockedView2x1( A02,
                       A12 );
    ABL.LockedView1x2( A20, A21 );
    ABR.LockedView( A22 );
#ifndef RELEASE
    PopCallStack();
#endif
}

#undef DM
#undef M

} // namespace elemental

#endif /* ELEMENTAL_PARTITIONING_HPP */

