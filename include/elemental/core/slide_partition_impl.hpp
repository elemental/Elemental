/*
   Copyright (c) 2009-2013, Jack Poulson
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
// SlidePartitionUp
//

template<typename T,typename Int>
inline void
SlidePartitionUp
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionUp [Matrix]");
#endif
    View( AT, A0 );
    View2x1( AB, A1,
                 A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
SlidePartitionUp
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionUp [DistMatrix]");
#endif
    View( AT, A0 );
    View2x1( AB, A1,
                 A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
SlideLockedPartitionUp
( M& AT, const M& A0,
         const M& A1,
  M& AB, const M& A2 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionUp [Matrix]");
#endif
    LockedView( AT, A0 );
    LockedView2x1( AB, A1,
                       A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
SlideLockedPartitionUp
( DM& AT, const DM& A0,
          const DM& A1,
  DM& AB, const DM& A2 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionUp [DistMatrix]");
#endif
    LockedView( AT, A0 );
    LockedView2x1( AB, A1,
                       A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// SlidePartitionDown
//

template<typename T,typename Int>
inline void
SlidePartitionDown
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionDown [Matrix]");
#endif
    View2x1( AT, A0,
                 A1 );
    View( AB, A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
SlidePartitionDown
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionDown [DistMatrix]");
#endif
    View2x1( AT, A0,
                 A1 );
    View( AB, A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
SlideLockedPartitionDown
( M& AT, const M& A0,
         const M& A1,
  M& AB, const M& A2 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionDown [Matrix]");
#endif
    LockedView2x1( AT, A0,
                       A1 );
    LockedView( AB, A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
SlideLockedPartitionDown
( DM& AT, const DM& A0,
          const DM& A1,
  DM& AB, const DM& A2 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionDown [DistMatrix]");
#endif
    LockedView2x1( AT, A0,
                       A1 );
    LockedView( AB, A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// SlidePartitionLeft
//

template<typename T,typename Int>
inline void
SlidePartitionLeft
( M& AL, M& AR,
  M& A0, M& A1, M& A2 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionLeft [Matrix]");
#endif
    View( AL, A0 );
    View1x2( AR, A1, A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
SlidePartitionLeft
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionLeft [DistMatrix]");
#endif
    View( AL, A0 );
    View1x2( AR, A1, A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
SlideLockedPartitionLeft
( M& AL, M& AR,
  const M& A0, const M& A1, const M& A2 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionLeft [Matrix]");
#endif
    LockedView( AL, A0 );
    LockedView1x2( AR, A1, A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
SlideLockedPartitionLeft
( DM& AL, DM& AR,
  const DM& A0, const DM& A1, const DM& A2 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionLeft [DistMatrix]");
#endif
    LockedView( AL, A0 );
    LockedView1x2( AR, A1, A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// SlidePartitionRight
//

template<typename T,typename Int>
inline void
SlidePartitionRight
( M& AL, M& AR,
  M& A0, M& A1, M& A2 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionRight [Matrix]");
#endif
    View1x2( AL, A0, A1 );
    View( AR, A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
SlidePartitionRight
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionRight [DistMatrix]");
#endif
    View1x2( AL, A0, A1 );
    View( AR, A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
SlideLockedPartitionRight
( M& AL, M& AR,
  const M& A0, const M& A1, const M& A2 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionRight [Matrix]");
#endif
    LockedView1x2( AL, A0, A1 );
    LockedView( AR, A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
SlideLockedPartitionRight
( DM& AL, DM& AR,
  const DM& A0, const DM& A1, const DM& A2 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionRight [DistMatrix]");
#endif
    LockedView1x2( AL, A0, A1 );
    LockedView( AR, A2 );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// SlidePartitionUpDiagonal
//

template<typename T,typename Int>
inline void
SlidePartitionUpDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionUpDiagonal [Matrix]");
#endif
    View( ATL, A00 );
    View1x2( ATR, A01, A02 );
    View2x1( ABL, A10,
                  A20 );
    View2x2( ABR, A11, A12,
                  A21, A22 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
SlidePartitionUpDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionUpDiagonal [DistMatrix]");
#endif
    View( ATL, A00 );
    View1x2( ATR, A01, A02 );
    View2x1( ABL, A10,
                  A20 );
    View2x2( ABR, A11, A12,
                  A21, A22 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
SlideLockedPartitionUpDiagonal
( M& ATL, M& ATR, const M& A00, const M& A01, const M& A02,
                  const M& A10, const M& A11, const M& A12,
  M& ABL, M& ABR, const M& A20, const M& A21, const M& A22 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionUpDiagonal [Matrix]");
#endif
    LockedView( ATL, A00 );
    LockedView1x2( ATR, A01, A02 );
    LockedView2x1( ABL, A10,
                        A20 );
    LockedView2x2( ABR, A11, A12,
                        A21, A22 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
SlideLockedPartitionUpDiagonal
( DM& ATL, DM& ATR, const DM& A00, const DM& A01, const DM& A02,
                    const DM& A10, const DM& A11, const DM& A12,
  DM& ABL, DM& ABR, const DM& A20, const DM& A21, const DM& A22 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionUpDiagonal [DistMatrix]");
#endif
    LockedView( ATL, A00 );
    LockedView1x2( ATR, A01, A02 );
    LockedView2x1( ABL, A10,
                        A20 );
    LockedView2x2( ABR, A11, A12,
                        A21, A22 );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// SlidePartitionDownDiagonal
//

template<typename T,typename Int>
inline void
SlidePartitionDownDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionDownDiagonal [Matrix]");
#endif
    View2x2( ATL, A00, A01,
                  A10, A11 );
    View2x1( ATR, A02,
                  A12 );
    View1x2( ABL, A20, A21 );
    View( ABR, A22 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
SlidePartitionDownDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22 )
{
#ifndef RELEASE
    PushCallStack("SlidePartitionDownDiagonal [DistMatrix]");
#endif
    View2x2( ATL, A00, A01,
                  A10, A11 );
    View2x1( ATR, A02,
                  A12 );
    View1x2( ABL, A20, A21 );
    View( ABR, A22 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
SlideLockedPartitionDownDiagonal
( M& ATL, M& ATR, const M& A00, const M& A01, const M& A02,
                  const M& A10, const M& A11, const M& A12,
  M& ABL, M& ABR, const M& A20, const M& A21, const M& A22 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionDownDiagonal [Matrix]");
#endif
    LockedView2x2( ATL, A00, A01,
                        A10, A11 );
    LockedView2x1( ATR, A02,
                        A12 );
    LockedView1x2( ABL, A20, A21 );
    LockedView( ABR, A22 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
SlideLockedPartitionDownDiagonal
( DM& ATL, DM& ATR, const DM& A00, const DM& A01, const DM& A02,
                    const DM& A10, const DM& A11, const DM& A12,
  DM& ABL, DM& ABR, const DM& A20, const DM& A21, const DM& A22 )
{
#ifndef RELEASE
    PushCallStack("SlideLockedPartitionDownDiagonal [DistMatrix]");
#endif
    LockedView2x2( ATL, A00, A01,
                        A10, A11 );
    LockedView2x1( ATR, A02,
                        A12 );
    LockedView1x2( ABL, A20, A21 );
    LockedView( ABR, A22 );
#ifndef RELEASE
    PopCallStack();
#endif
}

#undef DM
#undef M

} // namespace elem
