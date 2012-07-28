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
inline void
PartitionUp
( M& A, M& AT,
        M& AB, Int heightAB )
{
#ifndef RELEASE
    PushCallStack("PartitionUp [Matrix]");
    if( heightAB < 0 )
        throw std::logic_error
        ("Height of bottom partition must be non-negative");
#endif
    heightAB = std::min(heightAB,A.Height());
    const Int heightAT = A.Height()-heightAB;
    AT.View( A, 0,        0, heightAT, A.Width() );
    AB.View( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V,typename Int>
inline void
PartitionUp
( DM& A, DM& AT,
         DM& AB, Int heightAB )
{
#ifndef RELEASE
    PushCallStack("PartitionUp [DistMatrix]");
    if( heightAB < 0 )
        throw std::logic_error
        ("Height of bottom partition must be non-negative");
#endif
    heightAB = std::min(heightAB,A.Height());
    const Int heightAT = A.Height()-heightAB;
    AT.View( A, 0,        0, heightAT, A.Width() );
    AB.View( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
LockedPartitionUp
( const M& A, M& AT,
              M& AB, Int heightAB )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionUp [Matrix]");
    if( heightAB < 0 )
        throw std::logic_error
        ("Height of bottom partition must be non-negative");
#endif
    heightAB = std::min(heightAB,A.Height());
    const Int heightAT = A.Height()-heightAB;
    AT.LockedView( A, 0,        0, heightAT, A.Width() );
    AB.LockedView( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V,typename Int>
inline void
LockedPartitionUp
( const DM& A, DM& AT,
               DM& AB, Int heightAB )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionUp [DistMatrix]");
    if( heightAB < 0 )
        throw std::logic_error
        ("Height of bottom partition must be non-negative");
#endif
    heightAB = std::min(heightAB,A.Height());
    const Int heightAT = A.Height()-heightAB;
    AT.LockedView( A, 0,        0, heightAT, A.Width() );
    AB.LockedView( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// PartitionDown
//

template<typename T,typename Int>
inline void
PartitionDown
( M& A, M& AT,
        M& AB, Int heightAT ) 
{
#ifndef RELEASE
    PushCallStack("PartitionDown [Matrix]");
    if( heightAT < 0 )
        throw std::logic_error("Height of top partition must be non-negative");
#endif
    heightAT = std::min(heightAT,A.Height());
    const Int heightAB = A.Height()-heightAT;
    AT.View( A, 0,        0, heightAT, A.Width() );
    AB.View( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V,typename Int>
inline void
PartitionDown
( DM& A, DM& AT,
         DM& AB, Int heightAT )
{
#ifndef RELEASE
    PushCallStack("PartitionDown [DistMatrix]");
    if( heightAT < 0 )
        throw std::logic_error("Height of top partition must be non-negative");
#endif
    heightAT = std::min(heightAT,A.Height());
    const Int heightAB = A.Height()-heightAT;
    AT.View( A, 0,        0, heightAT, A.Width() );
    AB.View( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
LockedPartitionDown
( const M& A, M& AT,
              M& AB, Int heightAT ) 
{
#ifndef RELEASE
    PushCallStack("LockedPartitionDown [Matrix]");
    if( heightAT < 0 )
        throw std::logic_error("Height of top partition must be non-negative");
#endif
    heightAT = std::min(heightAT,A.Height());
    const Int heightAB = A.Height()-heightAT;
    AT.LockedView( A, 0,        0, heightAT, A.Width() );
    AB.LockedView( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V,typename Int>
inline void
LockedPartitionDown
( const DM& A, DM& AT,
               DM& AB, Int heightAT )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionDown [DistMatrix]");
    if( heightAT < 0 )
        throw std::logic_error("Height of top partition must be non-negative");
#endif
    heightAT = std::min(heightAT,A.Height());
    const Int heightAB = A.Height()-heightAT;
    AT.LockedView( A, 0,        0, heightAT, A.Width() );
    AB.LockedView( A, heightAT, 0, heightAB, A.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// PartitionLeft
//

template<typename T,typename Int>
inline void
PartitionLeft( M& A, M& AL, M& AR, Int widthAR )
{
#ifndef RELEASE
    PushCallStack("PartitionLeft [Matrix]");
    if( widthAR < 0 )
        throw std::logic_error("Width of right partition must be non-negative");
#endif
    widthAR = std::min(widthAR,A.Width());
    const Int widthAL = A.Width()-widthAR;
    AL.View( A, 0, 0,       A.Height(), widthAL );
    AR.View( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V,typename Int>
inline void
PartitionLeft( DM& A, DM& AL, DM& AR, Int widthAR )
{
#ifndef RELEASE
    PushCallStack("PartitionLeft [DistMatrix]");
    if( widthAR < 0 )
        throw std::logic_error("Width of right partition must be non-negative");
#endif
    widthAR = std::min(widthAR,A.Width());
    const Int widthAL = A.Width()-widthAR;
    AL.View( A, 0, 0,       A.Height(), widthAL );
    AR.View( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
LockedPartitionLeft( const M& A, M& AL, M& AR, Int widthAR )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionLeft [Matrix]");
    if( widthAR < 0 )
        throw std::logic_error("Width of right partition must be non-negative");
#endif
    widthAR = std::min(widthAR,A.Width());
    const Int widthAL = A.Width()-widthAR;
    AL.LockedView( A, 0, 0,       A.Height(), widthAL );
    AR.LockedView( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V,typename Int>
inline void
LockedPartitionLeft( const DM& A, DM& AL, DM& AR, Int widthAR )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionLeft [DistMatrix]");
    if( widthAR < 0 )
        throw std::logic_error("Width of right partition must be non-negative");
#endif
    widthAR = std::min(widthAR,A.Width());
    const Int widthAL = A.Width()-widthAR;
    AL.LockedView( A, 0, 0,       A.Height(), widthAL );
    AR.LockedView( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// PartitionRight
//

template<typename T,typename Int>
inline void
PartitionRight( M& A, M& AL, M& AR, Int widthAL )
{
#ifndef RELEASE
    PushCallStack("PartitionRight [Matrix]");
    if( widthAL < 0 )
        throw std::logic_error("Width of left partition must be non-negative");
#endif
    widthAL = std::min(widthAL,A.Width());
    const Int widthAR = A.Width()-widthAL;
    AL.View( A, 0, 0,       A.Height(), widthAL );
    AR.View( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V,typename Int>
inline void
PartitionRight( DM& A, DM& AL, DM& AR, Int widthAL )
{
#ifndef RELEASE
    PushCallStack("PartitionRight [DistMatrix]");
    if( widthAL < 0 )
        throw std::logic_error("Width of left partition must be non-negative");
#endif
    widthAL = std::min(widthAL,A.Width());
    const Int widthAR = A.Width()-widthAL;
    AL.View( A, 0, 0,       A.Height(), widthAL );
    AR.View( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
LockedPartitionRight( const M& A, M& AL, M& AR, Int widthAL )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionRight [Matrix]");
    if( widthAL < 0 )
        throw std::logic_error("Width of left partition must be non-negative");
#endif
    widthAL = std::min(widthAL,A.Width());
    const Int widthAR = A.Width()-widthAL;
    AL.LockedView( A, 0, 0,       A.Height(), widthAL );
    AR.LockedView( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V,typename Int>
inline void
LockedPartitionRight( const DM& A, DM& AL, DM& AR, Int widthAL )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionRight [DistMatrix]");
    if( widthAL < 0 )
        throw std::logic_error("Width of left partition must be non-negative");
#endif
    widthAL = std::min(widthAL,A.Width());
    const Int widthAR = A.Width()-widthAL;
    AL.LockedView( A, 0, 0,       A.Height(), widthAL );
    AR.LockedView( A, 0, widthAL, A.Height(), widthAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// PartitionUpDiagonal
//

template<typename T,typename Int>
inline void
PartitionUpDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagABR )
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

template<typename T, Distribution U, Distribution V,typename Int>
inline void
PartitionUpDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagABR )
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

template<typename T,typename Int>
inline void
LockedPartitionUpDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagABR )
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

template<typename T, Distribution U, Distribution V,typename Int>
inline void
LockedPartitionUpDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagABR )
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

//
// PartitionUpLeftDiagonal
//

template<typename T,typename Int>
inline void
PartitionUpLeftDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagABR )
{
#ifndef RELEASE
    PushCallStack("PartitionUpLeftDiagonal [Matrix]");
    if( diagABR < 0 )
        throw std::logic_error("Bottom-right size must be non-negative");
#endif
    const Int minDim = std::min( A.Height(), A.Width() );
    diagABR = std::min(diagABR,minDim);
    const Int sizeATL = minDim - diagABR;
    const Int remHeight = A.Height()-sizeATL;
    const Int remWidth = A.Width()-sizeATL;
    ATL.View( A, 0,       0,       sizeATL,   sizeATL  );
    ATR.View( A, 0,       sizeATL, sizeATL,   remWidth );
    ABL.View( A, sizeATL, 0,       remHeight, sizeATL  );
    ABR.View( A, sizeATL, sizeATL, remHeight, remWidth );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V,typename Int>
inline void
PartitionUpLeftDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagABR )
{
#ifndef RELEASE
    PushCallStack("PartitionUpLeftDiagonal [DistMatrix]");
    if( diagABR < 0 )
        throw std::logic_error("Bottom-right size must be non-negative");
#endif
    const Int minDim = std::min( A.Height(), A.Width() );
    diagABR = std::min(diagABR,minDim);
    const Int sizeATL = minDim - diagABR;
    const Int remHeight = A.Height()-sizeATL;
    const Int remWidth = A.Width()-sizeATL;
    ATL.View( A, 0,       0,       sizeATL,   sizeATL  );
    ATR.View( A, 0,       sizeATL, sizeATL,   remWidth );
    ABL.View( A, sizeATL, 0,       remHeight, sizeATL  );
    ABR.View( A, sizeATL, sizeATL, remHeight, remWidth );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
LockedPartitionUpLeftDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagABR )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionUpLeftDiagonal [Matrix]");
    if( diagABR < 0 )
        throw std::logic_error("Bottom-right size must be non-negative");
#endif
    const Int minDim = std::min( A.Height(), A.Width() );
    diagABR = std::min(diagABR,minDim);
    const Int sizeATL = minDim - diagABR;
    const Int remHeight = A.Height()-sizeATL;
    const Int remWidth = A.Width()-sizeATL;
    ATL.LockedView( A, 0,       0,       sizeATL,   sizeATL  );
    ATR.LockedView( A, 0,       sizeATL, sizeATL,   remWidth );
    ABL.LockedView( A, sizeATL, 0,       remHeight, sizeATL  );
    ABR.LockedView( A, sizeATL, sizeATL, remHeight, remWidth );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V,typename Int>
inline void
LockedPartitionUpLeftDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagABR )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionUpLeftDiagonal [DistMatrix]");
    if( diagABR < 0 )
        throw std::logic_error("Bottom-right size must be non-negative");
#endif
    const Int minDim = std::min( A.Height(), A.Width() );
    diagABR = std::min(diagABR,minDim);
    const Int sizeATL = minDim - diagABR;
    const Int remHeight = A.Height()-sizeATL;
    const Int remWidth = A.Width()-sizeATL;
    ATL.LockedView( A, 0,       0,       sizeATL,   sizeATL  );
    ATR.LockedView( A, 0,       sizeATL, sizeATL,   remWidth );
    ABL.LockedView( A, sizeATL, 0,       remHeight, sizeATL  );
    ABR.LockedView( A, sizeATL, sizeATL, remHeight, remWidth );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// PartitionUpRightDiagonal
//

template<typename T,typename Int>
inline void
PartitionUpRightDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagABR )
{
#ifndef RELEASE
    PushCallStack("PartitionUpRightDiagonal [Matrix]");
    if( diagABR < 0 )
        throw std::logic_error("Bottom-right size must be non-negative");
#endif
    const Int minDim = std::min(A.Height(),A.Width());
    diagABR = std::min(diagABR,minDim);
    const Int remHeight = A.Height()-diagABR;
    const Int remWidth = A.Width()-diagABR;
    ATL.View( A, 0,         0,        remHeight, remWidth );
    ATR.View( A, 0,         remWidth, remHeight, diagABR  );
    ABL.View( A, remHeight, 0,        diagABR,   remWidth );
    ABR.View( A, remHeight, remWidth, diagABR,   diagABR  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V,typename Int>
inline void
PartitionUpRightDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagABR )
{
#ifndef RELEASE
    PushCallStack("PartitionUpRightDiagonal [DistMatrix]");
    if( diagABR < 0 )
        throw std::logic_error("Bottom-right size must be non-negative");
#endif
    const Int minDim = std::min(A.Height(),A.Width());
    diagABR = std::min(diagABR,minDim);
    const Int remHeight = A.Height()-diagABR;
    const Int remWidth = A.Width()-diagABR;
    ATL.View( A, 0,         0,        remHeight, remWidth );
    ATR.View( A, 0,         remWidth, remHeight, diagABR  );
    ABL.View( A, remHeight, 0,        diagABR,   remWidth );
    ABR.View( A, remHeight, remWidth, diagABR,   diagABR  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
LockedPartitionUpRightDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagABR )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionUpRightDiagonal [Matrix]");
    if( diagABR < 0 )
        throw std::logic_error("Bottom-right size must be non-negative");
#endif
    const Int minDim = std::min(A.Height(),A.Width());
    diagABR = std::min(diagABR,minDim);
    const Int remHeight = A.Height()-diagABR;
    const Int remWidth = A.Width()-diagABR;
    ATL.View( A, 0,         0,        remHeight, remWidth );
    ATR.View( A, 0,         remWidth, remHeight, diagABR  );
    ABL.View( A, remHeight, 0,        diagABR,   remWidth );
    ABR.View( A, remHeight, remWidth, diagABR,   diagABR  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V,typename Int>
inline void
LockedPartitionUpRightDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagABR )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionUpRightDiagonal [DistMatrix]");
    if( diagABR < 0 )
        throw std::logic_error("Bottom-right size must be non-negative");
#endif
    const Int minDim = std::min(A.Height(),A.Width());
    diagABR = std::min(diagABR,minDim);
    const Int remHeight = A.Height()-diagABR;
    const Int remWidth = A.Width()-diagABR;
    ATL.View( A, 0,         0,        remHeight, remWidth );
    ATR.View( A, 0,         remWidth, remHeight, diagABR  );
    ABL.View( A, remHeight, 0,        diagABR,   remWidth );
    ABR.View( A, remHeight, remWidth, diagABR,   diagABR  );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// PartitionDownDiagonal
//

template<typename T,typename Int>
inline void
PartitionDownDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagATL )
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

template<typename T, Distribution U, Distribution V,typename Int>
inline void
PartitionDownDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagATL )
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

template<typename T,typename Int>
inline void
LockedPartitionDownDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagATL )
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

template<typename T, Distribution U, Distribution V,typename Int>
inline void
LockedPartitionDownDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagATL )
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

//
// PartitionDownLeftDiagonal
//

template<typename T,typename Int>
inline void
PartitionDownLeftDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagATL )
{
#ifndef RELEASE
    PushCallStack("PartitionDownLeftDiagonal [Matrix]");
    if( diagATL < 0 )
        throw std::logic_error("Top-left size must be non-negative");
#endif
    const Int minDim = std::min(A.Height(),A.Width());
    diagATL = std::min(diagATL,minDim);
    const Int heightABR = A.Height()-diagATL;
    const Int widthABR = A.Width()-diagATL;
    ATL.View( A, 0,       0,       diagATL,   diagATL  );
    ATR.View( A, 0,       diagATL, diagATL,   widthABR );
    ABL.View( A, diagATL, 0,       heightABR, diagATL  );
    ABR.View( A, diagATL, diagATL, heightABR, widthABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V,typename Int>
inline void
PartitionDownLeftDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagATL )
{
#ifndef RELEASE
    PushCallStack("PartitionDownLeftDiagonal [DistMatrix]");
    if( diagATL < 0 )
        throw std::logic_error("Top-left size must be non-negative");
#endif
    const Int minDim = std::min(A.Height(),A.Width());
    diagATL = std::min(diagATL,minDim);
    const Int heightABR = A.Height()-diagATL;
    const Int widthABR = A.Width()-diagATL;
    ATL.View( A, 0,       0,       diagATL,   diagATL  );
    ATR.View( A, 0,       diagATL, diagATL,   widthABR );
    ABL.View( A, diagATL, 0,       heightABR, diagATL  );
    ABR.View( A, diagATL, diagATL, heightABR, widthABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
LockedPartitionDownLeftDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagATL )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionDownLeftDiagonal [Matrix]");
    if( diagATL < 0 )
        throw std::logic_error("Top-left size must be non-negative");
#endif
    const Int minDim = std::min(A.Height(),A.Width());
    diagATL = std::min(diagATL,minDim);
    const Int heightABR = A.Height()-diagATL;
    const Int widthABR = A.Width()-diagATL;
    ATL.LockedView( A, 0,       0,       diagATL,   diagATL  );
    ATR.LockedView( A, 0,       diagATL, diagATL,   widthABR );
    ABL.LockedView( A, diagATL, 0,       heightABR, diagATL  );
    ABR.LockedView( A, diagATL, diagATL, heightABR, widthABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V,typename Int>
inline void
LockedPartitionDownLeftDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagATL )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionDownLeftDiagonal [DistMatrix]");
    if( diagATL < 0 )
        throw std::logic_error("Top-left size must be non-negative");
#endif
    const Int minDim = std::min(A.Height(),A.Width());
    diagATL = std::min(diagATL,minDim);
    const Int heightABR = A.Height()-diagATL;
    const Int widthABR = A.Width()-diagATL;
    ATL.LockedView( A, 0,       0,       diagATL,   diagATL  );
    ATR.LockedView( A, 0,       diagATL, diagATL,   widthABR );
    ABL.LockedView( A, diagATL, 0,       heightABR, diagATL  );
    ABR.LockedView( A, diagATL, diagATL, heightABR, widthABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// PartitionDownRightDiagonal
//

template<typename T,typename Int>
inline void
PartitionDownRightDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagATL )
{
#ifndef RELEASE
    PushCallStack("PartitionDownRightDiagonal [Matrix]");
    if( diagATL < 0 )
        throw std::logic_error("Top-left size must be non-negative");
#endif
    const Int minDim = std::min( A.Height(), A.Width() );
    diagATL = std::min(diagATL,minDim);
    const Int sizeABR = minDim-diagATL;
    const Int remHeight = A.Height()-sizeABR;
    const Int remWidth = A.Width()-sizeABR;
    ATL.View( A, 0,         0,        remHeight, remWidth );
    ATR.View( A, 0,         remWidth, remHeight, sizeABR  );
    ABL.View( A, remHeight, 0,        sizeABR,   remWidth );
    ABR.View( A, remHeight, remWidth, sizeABR,   sizeABR  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V,typename Int>
inline void
PartitionDownRightDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagATL )
{
#ifndef RELEASE
    PushCallStack("PartitionDownRightDiagonal [DistMatrix]");
    if( diagATL < 0 )
        throw std::logic_error("Top-left size must be non-negative");
#endif
    const Int minDim = std::min( A.Height(), A.Width() );
    diagATL = std::min(diagATL,minDim);
    const Int sizeABR = minDim-diagATL;
    const Int remHeight = A.Height()-sizeABR;
    const Int remWidth = A.Width()-sizeABR;
    ATL.View( A, 0,         0,        remHeight, remWidth );
    ATR.View( A, 0,         remWidth, remHeight, sizeABR  );
    ABL.View( A, remHeight, 0,        sizeABR,   remWidth );
    ABR.View( A, remHeight, remWidth, sizeABR,   sizeABR  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
LockedPartitionDownRightDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagATL )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionDownRightDiagonal [Matrix]");
    if( diagATL < 0 )
        throw std::logic_error("Top-left size must be non-negative");
#endif
    const Int minDim = std::min( A.Height(), A.Width() );
    diagATL = std::min(diagATL,minDim);
    const Int sizeABR = minDim-diagATL;
    const Int remHeight = A.Height()-sizeABR;
    const Int remWidth = A.Width()-sizeABR;
    ATL.View( A, 0,         0,        remHeight, remWidth );
    ATR.View( A, 0,         remWidth, remHeight, sizeABR  );
    ABL.View( A, remHeight, 0,        sizeABR,   remWidth );
    ABR.View( A, remHeight, remWidth, sizeABR,   sizeABR  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Distribution U, Distribution V,typename Int>
inline void
LockedPartitionDownRightDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagATL )
{
#ifndef RELEASE
    PushCallStack("LockedPartitionDownRightDiagonal [DistMatrix]");
    if( diagATL < 0 )
        throw std::logic_error("Top-left size must be non-negative");
#endif
    const Int minDim = std::min( A.Height(), A.Width() );
    diagATL = std::min(diagATL,minDim);
    const Int sizeABR = minDim-diagATL;
    const Int remHeight = A.Height()-sizeABR;
    const Int remWidth = A.Width()-sizeABR;
    ATL.View( A, 0,         0,        remHeight, remWidth );
    ATR.View( A, 0,         remWidth, remHeight, sizeABR  );
    ABL.View( A, remHeight, 0,        sizeABR,   remWidth );
    ABR.View( A, remHeight, remWidth, sizeABR,   sizeABR  );
#ifndef RELEASE
    PopCallStack();
#endif
}

#undef DM
#undef M

} // namespace elem
