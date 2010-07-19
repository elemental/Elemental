/*
   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Elemental.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ELEMENTAL_PARTITIONING_HPP
#define ELEMENTAL_PARTITIONING_HPP 1

#include "elemental/dist_matrix.hpp"

// To make our life easier. Undef'd at the bottom of the header
#define M  Matrix<T>
#define DM DistMatrix<T,U,V>

namespace elemental {

template<typename T>
void PartitionUp
( M& A, M& AT,
        M& AB, int heightAB );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void PartitionUp
( DM& A, DM& AT,
         DM& AB, int heightAB );

template<typename T>
void PartitionDown
( M& A, M& AT,
        M& AB, int heightAT );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void PartitionDown
( DM& A, DM& AT,
         DM& AB, int heightAT );

template<typename T>
void PartitionLeft
( M& A, M& AL, M& AR, int widthAR );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void PartitionLeft
( DM& A, DM& AL, DM& AR, int widthAR );

template<typename T>
void PartitionRight
( M& A, M& AL, M& AR, int widthAL );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void PartitionRight
( DM& A, DM& AL, DM& AR, int widthAL );

template<typename T>
void PartitionUpDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, int diagABR );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void PartitionUpDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, int diagABR );

template<typename T>
void PartitionDownDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, int diagATL );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void PartitionDownDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, int diagATL );

template<typename T>
void LockedPartitionUp
( const M& A, M& AT,
              M& AB, int heightAB );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void LockedPartitionUp
( const DM& A, DM& AT,
               DM& AB, int heightAB );

template<typename T>
void LockedPartitionDown
( const M& A, M& AT,
              M& AB, int heightAT );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void LockedPartitionDown
( const DM& A, DM& AT,
               DM& AB, int heightAT );

template<typename T>
void LockedPartitionLeft
( const M& A, M& AL, M& AR, int widthAR );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void LockedPartitionLeft
( const DM& A, DM& AL, DM& AR, int widthAR );

template<typename T>
void LockedPartitionRight
( const M& A, M& AL, M& AR, int widthAL );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void LockedPartitionRight
( const DM& A, DM& AL, DM& AR, int widthAL );

template<typename T>
void LockedPartitionUpDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, int diagABR );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void LockedPartitionUpDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, int diagABR );

template<typename T>
void LockedPartitionDownDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, int diagATL );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void LockedPartitionDownDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, int diagATL );

template<typename T>
void RepartitionUp
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void RepartitionUp
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2 );

template<typename T>
void RepartitionDown
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void RepartitionDown
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2 );

template<typename T>
void RepartitionLeft
( M& AL, M& AR,
  M& A0, M& A1, M& A2 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void RepartitionLeft
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2 );

template<typename T>
void RepartitionRight
( M& AL, M& AR,
  M& A0, M& A1, M& A2 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void RepartitionRight
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2 );

template<typename T>
void RepartitionUpDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void RepartitionUpDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22 );

template<typename T>
void RepartitionDownDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void RepartitionDownDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22 );

template<typename T>
void LockedRepartitionUp
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void LockedRepartitionUp
( const DM& AT, DM& A0,
                DM& A1,
  const DM& AB, DM& A2 );

template<typename T>
void LockedRepartitionDown
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void LockedRepartitionDown
( const DM& AT, DM& A0,
                DM& A1,
  const DM& AB, DM& A2 );

template<typename T>
void LockedRepartitionLeft
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2      );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void LockedRepartitionLeft
( const DM& AL, const DM& AR,
  DM& A0, DM& A1, DM& A2     );

template<typename T>
void LockedRepartitionRight
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2      );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void LockedRepartitionRight
( const DM& AL, const DM& AR,
  DM& A0, DM& A1, DM& A2     );

template<typename T>
void LockedRepartitionUpDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void LockedRepartitionUpDiagonal
( const DM& ATL, const DM& ATR, DM& A00, DM& A01, DM& A02,
                                DM& A10, DM& A11, DM& A12,
  const DM& ABL, const DM& ABR, DM& A20, DM& A21, DM& A22 );

template<typename T>
void LockedRepartitionDownDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void LockedRepartitionDownDiagonal
( const DM& ATL, const DM& ATR, DM& A00, DM& A01, DM& A02,
                                DM& A10, DM& A11, DM& A12,
  const DM& ABL, const DM& ABR, DM& A20, DM& A21, DM& A22 );

template<typename T>
void SlidePartitionUp
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void SlidePartitionUp
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2 );

template<typename T>
void SlidePartitionDown
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void SlidePartitionDown
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2 );

template<typename T>
void SlidePartitionLeft
( M& AL, M& AR,
  M& A0, M& A1, M& A2 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void SlidePartitionLeft
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2 );

template<typename T>
void SlidePartitionRight
( M& AL, M& AR,
  M& A0, M& A1, M& A2 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void SlidePartitionRight
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2 );

template<typename T>
void SlidePartitionUpDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void SlidePartitionUpDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22 );

template<typename T>
void SlidePartitionDownDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void SlidePartitionDownDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22 );

template<typename T>
void SlideLockedPartitionUp
( M& AT, const M& A0,
         const M& A1,
  M& AB, const M& A2 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void SlideLockedPartitionUp
( DM& AT, const DM& A0,
          const DM& A1,
  DM& AB, const DM& A2 );

template<typename T>
void SlideLockedPartitionDown
( M& AT, const M& A0,
         const M& A1,
  M& AB, const M& A2 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void SlideLockedPartitionDown
( DM& AT, const DM& A0,
          const DM& A1,
  DM& AB, const DM& A2 );

template<typename T>
void SlideLockedPartitionLeft
( M& AL, M& AR,
  const M& A0, const M& A1, const M& A2 ); 

template<typename T, elemental::Distribution U, elemental::Distribution V>
void SlideLockedPartitionLeft
( DM& AL, DM& AR,
  const DM& A0, const DM& A1, const DM& A2 );

template<typename T>
void SlideLockedPartitionRight
( M& AL, M& AR,
  const M& A0, const M& A1, const M& A2 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void SlideLockedPartitionRight
( DM& AL, DM& AR,
  const DM& A0, const DM& A1, const DM& A2 );

template<typename T>
void SlideLockedPartitionUpDiagonal
( M& ATL, M& ATR, const M& A00, const M& A01, const M& A02,
                  const M& A10, const M& A11, const M& A12,
  M& ABL, M& ABR, const M& A20, const M& A21, const M& A22 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void SlideLockedPartitionUpDiagonal
( DM& ATL, DM& ATR, const DM& A00, const DM& A01, const DM& A02,
                    const DM& A10, const DM& A11, const DM& A12,
  DM& ABL, DM& ABR, const DM& A20, const DM& A21, const DM& A22 );

template<typename T>
void SlideLockedPartitionDownDiagonal
( M& ATL, M& ATR, const M& A00, const M& A01, const M& A02,
                  const M& A10, const M& A11, const M& A12,
  M& ABL, M& ABR, const M& A20, const M& A21, const M& A22 );

template<typename T, elemental::Distribution U, elemental::Distribution V>
void SlideLockedPartitionDownDiagonal
( DM& ATL, DM& ATR, const DM& A00, const DM& A01, const DM& A02,
                    const DM& A10, const DM& A11, const DM& A12,
  DM& ABL, DM& ABR, const DM& A20, const DM& A21, const DM& A22 );

} // elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename T>
inline void
elemental::PartitionUp
( M& A, M& AT,
        M& AB, int heightAB )
{
#ifndef RELEASE
    PushCallStack("PartitionUp [Matrix]");
    if( heightAB < 0 )
        throw std::logic_error
        ( "Height of bottom partition must be non-negative." );
    if( heightAB > A.Height() )
        throw std::logic_error( "Height of bottom partition is too large." );
#endif
    const int heightAT = A.Height()-heightAB;
    AT.View( A, 0,        0, heightAT, A.Width() );
    AB.View( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::PartitionUp
( DM& A, DM& AT,
         DM& AB, int heightAB )
{
#ifndef RELEASE
    PushCallStack("PartitionUp [DistMatrix]");
    if( heightAB < 0 )
        throw std::logic_error
        ( "Height of bottom partition must be non-negative." );
    if( heightAB > A.Height() )
        throw std::logic_error( "Height of bottom partition is too large." );
#endif
    const int heightAT = A.Height()-heightAB;
    AT.View( A, 0,        0, heightAT, A.Width() );
    AB.View( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::PartitionDown
( M& A, M& AT,
        M& AB, int heightAT ) 
{
#ifndef RELEASE
    PushCallStack("PartitionDown [Matrix]");
    if( heightAT < 0 )
        throw std::logic_error
        ( "Height of top partition must be non-negative." );
    if( heightAT > A.Height() )
        throw std::logic_error( "Height of top partition is too large." );
#endif
    const int heightAB = A.Height()-heightAT;
    AT.View( A, 0,        0, heightAT, A.Width() );
    AB.View( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::PartitionDown
( DM& A, DM& AT,
         DM& AB, int heightAT )
{
#ifndef RELEASE
    PushCallStack("PartitionDown [DistMatrix]");
    if( heightAT < 0 )
        throw std::logic_error
        ( "Height of top partition must be non-negative." );
    if( heightAT > A.Height() )
        throw std::logic_error( "Height of top partition is too large." );
#endif
    const int heightAB = A.Height()-heightAT;
    AT.View( A, 0,        0, heightAT, A.Width() );
    AB.View( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::PartitionLeft
( M& A, M& AL, M& AR, int widthAR )
{
#ifndef RELEASE
    PushCallStack("PartitionLeft [Matrix]");
    if( widthAR < 0 )
        throw std::logic_error
        ( "Width of right partition must be non-negative." );
    if( widthAR > A.Width() )
        throw std::logic_error( "Width of right partition is too large." );
#endif
    const int widthAL = A.Width()-widthAR;
    AL.View( A, 0, 0,       A.Height(), widthAL );
    AR.View( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::PartitionLeft
( DM& A, DM& AL, DM& AR, int widthAR )
{
#ifndef RELEASE
    PushCallStack("PartitionLeft [DistMatrix]");
    if( widthAR < 0 )
        throw std::logic_error
        ( "Width of right partition must be non-negative." );
    if( widthAR > A.Width() )
        throw std::logic_error( "Width of right partition is too large." );
#endif
    const int widthAL = A.Width()-widthAR;
    AL.View( A, 0, 0,       A.Height(), widthAL );
    AR.View( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::PartitionRight
( M& A, M& AL, M& AR, int widthAL )
{
#ifndef RELEASE
    PushCallStack("PartitionRight [Matrix]");
    if( widthAL < 0 )
        throw std::logic_error
        ( "Width of left partition must be non-negative." );
    if( widthAL > A.Width() )
        throw std::logic_error( "Width of left partition is too large." );
#endif
    const int widthAR = A.Width()-widthAL;
    AL.View( A, 0, 0,       A.Height(), widthAL );
    AR.View( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::PartitionRight
( DM& A, DM& AL, DM& AR, int widthAL )
{
#ifndef RELEASE
    PushCallStack("PartitionRight [DistMatrix]");
    if( widthAL < 0 )
        throw std::logic_error
        ( "Width of left partition must be non-negative." );
    if( widthAL > A.Width() )
        throw std::logic_error( "Width of left partition is too large." );
#endif
    const int widthAR = A.Width()-widthAL;
    AL.View( A, 0, 0,       A.Height(), widthAL );
    AR.View( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::PartitionUpDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, int diagABR )
{
#ifndef RELEASE
    PushCallStack("PartitionUpDiagonal [Matrix]");
    if( diagABR < 0 )
        throw std::logic_error( "Bottom-right size must be non-negative." );
    if( diagABR > A.Height() || diagABR > A.Width() )
        throw std::logic_error( "Bottom-right size is too large." );
#endif
    const int minDim = std::min( A.Height(), A.Width() );
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

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::PartitionUpDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, int diagABR )
{
#ifndef RELEASE
    PushCallStack("PartitionUpDiagonal [DistMatrix]");
    if( diagABR < 0 )
        throw std::logic_error( "Bottom-right size must be non-negative." );
    if( diagABR > A.Height() || diagABR > A.Width() )
        throw std::logic_error( "Bottom-right size is too large." );
#endif
    const int minDim = std::min( A.Height(), A.Width() );
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
elemental::PartitionDownDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, int diagATL )
{
#ifndef RELEASE
    PushCallStack("PartitionDownDiagonal [Matrix]");
    if( diagATL < 0 )
        throw std::logic_error( "Top-left size must be non-negative." );
    if( diagATL > A.Height() || diagATL > A.Width() )
        throw std::logic_error( "Top-left size is too large." );
#endif
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

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::PartitionDownDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, int diagATL )
{
#ifndef RELEASE
    PushCallStack("PartitionDownDiagonal [DistMatrix]");
    if( diagATL < 0 )
        throw std::logic_error( "Top-left size must be non-negative." );
    if( diagATL > A.Height() || diagATL > A.Width() )
        throw std::logic_error( "Top-left size is too large." );
#endif
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
elemental::LockedPartitionUp
( const M& A, M& AT,
              M& AB, int heightAB )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionUp [Matrix]");
    if( heightAB < 0 )
        throw std::logic_error
        ( "Height of bottom partition must be non-negative." );
    if( heightAB > A.Height() )
        throw std::logic_error( "Height of bottom partition is too large." );
#endif
    const int heightAT = A.Height()-heightAB;
    AT.LockedView( A, 0,        0, heightAT, A.Width() );
    AB.LockedView( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::LockedPartitionUp
( const DM& A, DM& AT,
               DM& AB, int heightAB )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionUp [DistMatrix]");
    if( heightAB < 0 )
        throw std::logic_error
        ( "Height of bottom partition must be non-negative." );
    if( heightAB > A.Height() )
        throw std::logic_error( "height of bottom partition is too large." );
#endif
    const int heightAT = A.Height()-heightAB;
    AT.LockedView( A, 0,        0, heightAT, A.Width() );
    AB.LockedView( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::LockedPartitionDown
( const M& A, M& AT,
              M& AB, int heightAT ) 
{
#ifndef RELEASE
    PushCallStack("LockedPartitionDown [Matrix]");
    if( heightAT < 0 )
        throw std::logic_error
        ( "Height of top partition must be non-negative." );
    if( heightAT > A.Height() )
        throw std::logic_error( "Height of bottom partition is too large." );
#endif
    const int heightAB = A.Height()-heightAT;
    AT.LockedView( A, 0,        0, heightAT, A.Width() );
    AB.LockedView( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::LockedPartitionDown
( const DM& A, DM& AT,
               DM& AB, int heightAT )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionDown [DistMatrix]");
    if( heightAT < 0 )
        throw std::logic_error
        ( "Height of top partition must be non-negative." );
    if( heightAT > A.Height() )
        throw std::logic_error( "Height of top partition is too large." );
#endif
    const int heightAB = A.Height()-heightAT;
    AT.LockedView( A, 0,        0, heightAT, A.Width() );
    AB.LockedView( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::LockedPartitionLeft
( const M& A, M& AL, M& AR, int widthAR )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionLeft [Matrix]");
    if( widthAR < 0 )
        throw std::logic_error
        ( "Width of right partition must be non-negative." );
    if( widthAR > A.Width() )
        throw std::logic_error( "Width of right partition is too large." );
#endif
    const int widthAL = A.Width()-widthAR;
    AL.LockedView( A, 0, 0,       A.Height(), widthAL );
    AR.LockedView( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::LockedPartitionLeft
( const DM& A, DM& AL, DM& AR, int widthAR )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionLeft [DistMatrix]");
    if( widthAR < 0 )
        throw std::logic_error
        ( "Width of right partition must be non-negative." );
    if( widthAR > A.Width() )
        throw std::logic_error( "Width of right partition is too large." );
#endif
    const int widthAL = A.Width()-widthAR;
    AL.LockedView( A, 0, 0,       A.Height(), widthAL );
    AR.LockedView( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::LockedPartitionRight
( const M& A, M& AL, M& AR, int widthAL )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionRight [Matrix]");
    if( widthAL < 0 )
        throw std::logic_error
        ( "Width of left partition must be non-negative." );
    if( widthAL > A.Width() )
        throw std::logic_error( "Width of left partition is too large." );
#endif
    const int widthAR = A.Width()-widthAL;
    AL.LockedView( A, 0, 0,       A.Height(), widthAL );
    AR.LockedView( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::LockedPartitionRight
( const DM& A, DM& AL, DM& AR, int widthAL )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionRight [DistMatrix]");
    if( widthAL < 0 )
        throw std::logic_error
        ( "Width of left partition must be non-negative." );
    if( widthAL > A.Width() )
        throw std::logic_error( "Width of left partition is too large." );
#endif
    const int widthAR = A.Width()-widthAL;
    AL.LockedView( A, 0, 0,       A.Height(), widthAL );
    AR.LockedView( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::LockedPartitionUpDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, int diagABR )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionUpDiagonal [Matrix]");
    if( diagABR < 0 )
        throw std::logic_error( "Bottom-right size must be non-negative." );
    if( diagABR > A.Height() || diagABR > A.Width() )
        throw std::logic_error( "Bottom-right size is too large." );
#endif
    const int minDim = std::min( A.Height(), A.Width() );
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

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::LockedPartitionUpDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, int diagABR )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionUpDiagonal [DistMatrix]");
    if( diagABR < 0 )
        throw std::logic_error( "Bottom-right size must be non-negative." );
    if( diagABR > A.Height() || diagABR > A.Width() )
        throw std::logic_error( "Bottom-right size is too large." );
#endif
    const int minDim = std::min( A.Height(), A.Width() );
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
elemental::LockedPartitionDownDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, int diagATL )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionDownDiagonal [Matrix]");
    if( diagATL < 0 )
        throw std::logic_error( "Top-left size must be non-negative." );
    if( diagATL > A.Height() || diagATL > A.Width() )
        throw std::logic_error( "Top-left size is too large." );
#endif
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

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::LockedPartitionDownDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, int diagATL )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionDownDiagonal [DistMatrix]");
    if( diagATL < 0 )
        throw std::logic_error( "Top-left size must be non-negative." );
    if( diagATL > A.Height() || diagATL > A.Width() )
        throw std::logic_error( "Top-left size is too large." );
#endif
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
elemental::RepartitionUp
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2 )
{
#ifndef RELEASE
    PushCallStack("RepartitionUp [Matrix]");
    if( (AT.Buffer() + AT.Height()) != AB.Buffer() )
        throw std::logic_error( "Noncontiguous 2x1 array of matrices." );
#endif
    int bsize = std::min( AT.Height(), Blocksize() );
    int offset = AT.Height()-bsize; 
    A0.View( AT, 0,      0, offset, AT.Width() );
    A1.View( AT, offset, 0, bsize,  AT.Width() );
    A2.View( AB );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::RepartitionUp
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2 )
{
#ifndef RELEASE
    PushCallStack("RepartitionUp [DistMatrix]");
    if( (AT.LocalMatrix().Buffer() + AT.LocalHeight()) != 
         AB.LocalMatrix().Buffer()                        )
    {
        throw std::logic_error
        ( "Noncontiguous 2x1 array of distributed matrices." );
    }
#endif
    int bsize = std::min( AT.Height(), Blocksize() );
    int offset = AT.Height()-bsize; 
    A0.View( AT, 0,      0, offset, AT.Width() );
    A1.View( AT, offset, 0, bsize,  AT.Width() );
    A2.View( AB );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::RepartitionDown
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2 )
{
#ifndef RELEASE
    PushCallStack("RepartitionDown [Matrix]");
    if( (AT.Buffer() + AT.Height()) != AB.Buffer() )
        throw std::logic_error( "Noncontiguous 2x1 array of matrices." );
#endif
    int bsize = std::min( AB.Height(), Blocksize() );
    int offset = AB.Height()-bsize; 
    A0.View( AT );
    A1.View( AB, 0,     0, bsize,  AB.Width() );
    A2.View( AB, bsize, 0, offset, AB.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::RepartitionDown
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2 )
{
#ifndef RELEASE
    PushCallStack("RepartitionDown [DistMatrix]");
    if( (AT.LocalMatrix().Buffer() + AT.LocalHeight()) != 
         AB.LocalMatrix().Buffer() )
    {
        throw std::logic_error
        ( "Noncontiguous 2x1 array of distributed matrices." );
    }
#endif
    int bsize = std::min( AB.Height(), Blocksize() );
    int offset = AB.Height()-bsize; 
    A0.View( AT );
    A1.View( AB, 0,     0, bsize,  AB.Width() );
    A2.View( AB, bsize, 0, offset, AB.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::RepartitionLeft
( M& AL, M& AR,
  M& A0, M& A1, M& A2 )
{
#ifndef RELEASE
    PushCallStack("RepartitionLeft [Matrix]");
    if( (AL.Buffer() + AL.Width()*AL.LDim()) != AR.Buffer() )
    {
        throw std::logic_error( "Noncontiguous 1x2 array of matrices." );
    }
#endif
    int bsize = std::min( AL.Width(), Blocksize() );
    int offset = AL.Width()-bsize;
    A0.View( AL, 0, 0,      AL.Height(), offset );
    A1.View( AL, 0, offset, AL.Height(), bsize  );
    A2.View( AR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::RepartitionLeft
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2 )
{
#ifndef RELEASE
    PushCallStack("RepartitionLeft [DistMatrix]");
    if( (AL.LocalMatrix().Buffer() + AL.LocalWidth()*AL.LocalLDim())
         != AR.LocalMatrix().Buffer() )
    {
        throw std::logic_error
        ( "Noncontiguous 1x2 array of distributed matrices." );
    }
#endif
    int bsize = std::min( AL.Width(), Blocksize() );
    int offset = AL.Width()-bsize;
    A0.View( AL, 0, 0,      AL.Height(), offset );
    A1.View( AL, 0, offset, AL.Height(), bsize  );
    A2.View( AR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::RepartitionRight
( M& AL, M& AR,
  M& A0, M& A1, M& A2 )
{
#ifndef RELEASE
    PushCallStack("RepartitionRight [Matrix]");
    if( (AL.Buffer() + AL.Width()*AL.LDim()) != AR.Buffer() )
    {
        throw std::logic_error( "Noncontiguous 1x2 array of matrices." );
    }
#endif
    int bsize = std::min( AR.Width(), Blocksize() );
    int offset = AR.Width()-bsize;
    A0.View( AL );
    A1.View( AR, 0, 0,     AR.Height(), bsize  );
    A2.View( AR, 0, bsize, AR.Height(), offset );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::RepartitionRight
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2 )
{
#ifndef RELEASE
    PushCallStack("RepartitionRight [DistMatrix]");
    if( (AL.LocalMatrix().Buffer() + AL.LocalWidth()*AL.LocalLDim()) 
         != AR.LocalMatrix().Buffer() )
    {
        throw std::logic_error
        ( "Noncontiguous 1x2 array of distributed matrices." );
    }
#endif
    int bsize = std::min( AR.Width(), Blocksize() );
    int offset = AR.Width()-bsize;
    A0.View( AL );
    A1.View( AR, 0, 0,     AR.Height(), bsize  );
    A2.View( AR, 0, bsize, AR.Height(), offset );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::RepartitionUpDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22 )
{
#ifndef RELEASE
    PushCallStack("RepartitionUpDiagonal [Matrix]");
    if( (ATL.Buffer() + ATL.Height()) != ABL.Buffer() ||
        (ATR.Buffer() + ATR.Height()) != ABR.Buffer() ||
        (ATL.Buffer() + ATL.Width()*ATL.LDim()) != ATR.Buffer() ||
        (ABL.Buffer() + ABL.Width()*ABL.LDim()) != ABR.Buffer()    )
        throw std::logic_error( "Noncontiguous 2x2 grid of matrices." );
#endif
    int bsize = std::min( Blocksize(), std::min( ATL.Height(),
                                                 ATL.Width() ) );
    int vOffset = ATL.Height()-bsize;
    int hOffset = ATL.Width()-bsize;
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

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::RepartitionUpDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22 )
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
    {
        throw std::logic_error
        ( "Noncontiguous 2x2 grid of distributed matrices." );
    }
#endif
    int bsize = std::min( Blocksize(), std::min( ATL.Height(),
                                                 ATL.Width() ) );
    int vOffset = ATL.Height()-bsize;
    int hOffset = ATL.Width()-bsize;
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
elemental::RepartitionDownDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22 )
{
#ifndef RELEASE
    PushCallStack("RepartitionDownDiagonal [Matrix]");
    if( (ATL.Buffer() + ATL.Height()) != ABL.Buffer() ||
        (ATR.Buffer() + ATR.Height()) != ABR.Buffer() ||
        (ATL.Buffer() + ATL.Width()*ATL.LDim()) != ATR.Buffer() ||
        (ABL.Buffer() + ABL.Width()*ABL.LDim()) != ABR.Buffer()    )
    {
        throw std::logic_error( "Noncontiguous 2x2 grid of matrices." );
    }
#endif
    int bsize = std::min( Blocksize(), std::min( ABR.Height(),
                                                 ABR.Width() ) );
    int vOffset = ABR.Height()-bsize;
    int hOffset = ABR.Width()-bsize;
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

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::RepartitionDownDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22 )
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
    {
        throw std::logic_error
        ( "Noncontiguous 2x2 grid of distributed matrices." );
    }
#endif
    int bsize = std::min( Blocksize(), std::min( ABR.Height(),
                                                 ABR.Width() ) );
    int vOffset = ABR.Height()-bsize;
    int hOffset = ABR.Width()-bsize;
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
elemental::LockedRepartitionUp
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2 )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionUp [Matrix]");
    if( (AT.LockedBuffer() + AT.Height()) != AB.LockedBuffer() )
        throw std::logic_error( "Noncontiguous 2x1 array of matrices." );
#endif
    int bsize = std::min( AT.Height(), Blocksize() );
    int offset = AT.Height()-bsize;
    A0.LockedView( AT, 0,      0, offset, AT.Width() );
    A1.LockedView( AT, offset, 0, bsize,  AT.Width() );
    A2.LockedView( AB );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::LockedRepartitionUp
( const DM& AT, DM& A0,
                DM& A1,
  const DM& AB, DM& A2 )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionUp [DistMatrix]");
    if( (AT.LockedLocalMatrix().LockedBuffer() + AT.LocalHeight()) != 
         AB.LockedLocalMatrix().LockedBuffer() )
    {
        throw std::logic_error
        ( "Noncontiguous 2x1 array of distributed matrices." );
    }
#endif
    int bsize = std::min( AT.Height(), Blocksize() );
    int offset = AT.Height()-bsize;
    A0.LockedView( AT, 0,      0, offset, AT.Width() );
    A1.LockedView( AT, offset, 0, bsize,  AT.Width() );
    A2.LockedView( AB );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::LockedRepartitionDown
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2 )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionDown [Matrix]");
    if( (AT.LockedBuffer() + AT.Height()) != AB.LockedBuffer() )
        throw std::logic_error( "Noncontiguous 2x1 array of matrices." );
#endif
    int bsize = std::min( AB.Height(), Blocksize() );
    int offset = AB.Height()-bsize;
    A0.LockedView( AT );
    A1.LockedView( AB, 0,     0, bsize,  AB.Width() );
    A2.LockedView( AB, bsize, 0, offset, AB.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::LockedRepartitionDown
( const DM& AT, DM& A0,
                DM& A1,
  const DM& AB, DM& A2 )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionDown [DistMatrix]");
    if( (AT.LockedLocalMatrix().LockedBuffer() + AT.LocalHeight()) != 
         AB.LockedLocalMatrix().LockedBuffer() )
    {
        throw std::logic_error
        ( "Noncontiguous 2x1 array of distributed matrices." );
    }
#endif
    int bsize = std::min( AB.Height(), Blocksize() );
    int offset = AB.Height()-bsize;
    A0.LockedView( AT );
    A1.LockedView( AB, 0,     0, bsize,  AB.Width() );
    A2.LockedView( AB, bsize, 0, offset, AB.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::LockedRepartitionLeft
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2      )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionLeft [Matrix]");
    if( (AL.LockedBuffer() + AL.Width()*AL.LDim()) != AR.LockedBuffer() )
    {
        throw std::logic_error( "Noncontiguous 1x2 array of matrices." );
    }
#endif
    int bsize = std::min( AL.Width(), Blocksize() );
    int offset = AL.Width()-bsize;
    A0.LockedView( AL, 0, 0,      AL.Height(), offset );
    A1.LockedView( AL, 0, offset, AL.Height(), bsize  );
    A2.LockedView( AR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::LockedRepartitionLeft
( const DM& AL, const DM& AR,
  DM& A0, DM& A1, DM& A2     )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionLeft [DistMatrix]");
    if( (AL.LockedLocalMatrix().LockedBuffer() + AL.LocalWidth()*AL.LocalLDim())
         != AR.LockedLocalMatrix().LockedBuffer() )
    {
        throw std::logic_error
        ( "Noncontiguous 1x1 array of distributed matrices." );
    }
#endif
    int bsize = std::min( AL.Width(), Blocksize() );
    int offset = AL.Width()-bsize;
    A0.LockedView( AL, 0, 0,      AL.Height(), offset );
    A1.LockedView( AL, 0, offset, AL.Height(), bsize  );
    A2.LockedView( AR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::LockedRepartitionRight
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2      )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionRight [Matrix]");
    if( (AL.LockedBuffer() + AL.Width()*AL.LDim()) != AR.LockedBuffer() )
    {
        throw std::logic_error( "Noncontiguous 1x2 array of matrices." );
    }
#endif
    int bsize = std::min( AR.Width(), Blocksize() );
    int offset = AR.Width()-bsize;
    A0.LockedView( AL );
    A1.LockedView( AR, 0, 0,     AR.Height(), bsize  );
    A2.LockedView( AR, 0, bsize, AR.Height(), offset );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::LockedRepartitionRight
( const DM& AL, const DM& AR,
  DM& A0, DM& A1, DM& A2     )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionRight [DistMatrix]");
    if( (AL.LockedLocalMatrix().LockedBuffer() + AL.LocalWidth()*AL.LocalLDim())
         != AR.LockedLocalMatrix().LockedBuffer() )
    {
        throw std::logic_error
        ( "Noncontiguous 1x2 DistMatrices in LockedRepartitionRight." );
    }
#endif
    int bsize = std::min( AR.Width(), Blocksize() );
    int offset = AR.Width()-bsize;
    A0.LockedView( AL );
    A1.LockedView( AR, 0, 0,     AR.Height(), bsize  );
    A2.LockedView( AR, 0, bsize, AR.Height(), offset );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::LockedRepartitionUpDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22 )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionUpDiagonal [Matrix]");
    if( (ATL.LockedBuffer() + ATL.Height()) != ABL.LockedBuffer() ||
        (ATR.LockedBuffer() + ATR.Height()) != ABR.LockedBuffer() ||
        (ATL.LockedBuffer() + ATL.Width()*ATL.LDim()) != ATR.LockedBuffer() ||
        (ABL.LockedBuffer() + ABL.Width()*ABL.LDim()) != ABR.LockedBuffer()    )
    {
        throw std::logic_error( "Noncontiguous 2x2 grid of matrices." );
    }
#endif
    int bsize = std::min( Blocksize(), std::min( ATL.Height(),
                                                 ATL.Width() ) );
    int vOffset = ATL.Height()-bsize;
    int hOffset = ATL.Width()-bsize;
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

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::LockedRepartitionUpDiagonal
( const DM& ATL, const DM& ATR, DM& A00, DM& A01, DM& A02,
                                DM& A10, DM& A11, DM& A12,
  const DM& ABL, const DM& ABR, DM& A20, DM& A21, DM& A22 )
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
    {
        throw std::logic_error
        ( "Noncontiguous 2x2 grid of distributed matrices." );
    }
#endif
    int bsize = std::min( Blocksize(), std::min( ATL.Height(),
                                                 ATL.Width() ) );
    int vOffset = ATL.Height()-bsize;
    int hOffset = ATL.Width()-bsize;
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
elemental::LockedRepartitionDownDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22 )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionDownDiagonal [Matrix]");
    if( (ATL.LockedBuffer() + ATL.Height()) != ABL.LockedBuffer() ||
        (ATR.LockedBuffer() + ATR.Height()) != ABR.LockedBuffer() ||
        (ATL.LockedBuffer() + ATL.Width()*ATL.LDim()) != ATR.LockedBuffer() ||
        (ABL.LockedBuffer() + ABL.Width()*ABL.LDim()) != ABR.LockedBuffer()    )
    {
        throw std::logic_error( "Noncontiguous 2x2 grid of matrices." );
    }
#endif
    int bsize = std::min( Blocksize(), std::min( ABR.Height(),
                                                 ABR.Width() ) );
    int vOffset = ABR.Height()-bsize;
    int hOffset = ABR.Width()-bsize;
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

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::LockedRepartitionDownDiagonal
( const DM& ATL, const DM& ATR, DM& A00, DM& A01, DM& A02,
                                DM& A10, DM& A11, DM& A12,
  const DM& ABL, const DM& ABR, DM& A20, DM& A21, DM& A22 )
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
    {
        throw std::logic_error
        ( "Noncontiguous 2x2 grid of distributed matrices." );
    }
#endif
    int bsize = std::min( Blocksize(), std::min( ABR.Height(),
                                                 ABR.Width() ) );
    int vOffset = ABR.Height()-bsize;
    int hOffset = ABR.Width()-bsize;
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
elemental::SlidePartitionUp
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

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::SlidePartitionUp
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
elemental::SlidePartitionDown
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

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::SlidePartitionDown
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
elemental::SlidePartitionLeft
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

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::SlidePartitionLeft
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
elemental::SlidePartitionRight
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

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::SlidePartitionRight
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
elemental::SlidePartitionUpDiagonal
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

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::SlidePartitionUpDiagonal
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
elemental::SlidePartitionDownDiagonal
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

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::SlidePartitionDownDiagonal
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
elemental::SlideLockedPartitionUp
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

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::SlideLockedPartitionUp
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
elemental::SlideLockedPartitionDown
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

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::SlideLockedPartitionDown
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
elemental::SlideLockedPartitionLeft
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

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::SlideLockedPartitionLeft
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
elemental::SlideLockedPartitionRight
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

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::SlideLockedPartitionRight
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
elemental::SlideLockedPartitionUpDiagonal
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

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::SlideLockedPartitionUpDiagonal
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
elemental::SlideLockedPartitionDownDiagonal
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

template<typename T, elemental::Distribution U, elemental::Distribution V>
inline void
elemental::SlideLockedPartitionDownDiagonal
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

#endif /* ELEMENTAL_PARTITIONING_HPP */

