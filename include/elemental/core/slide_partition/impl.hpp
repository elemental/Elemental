/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_SLIDEPARTITION_IMPL_HPP
#define ELEM_CORE_SLIDEPARTITION_IMPL_HPP

namespace elem {

// To make our life easier. Undef'd at the bottom of the header
#define M  Matrix<T>
#define DM DistMatrix<T,U,V>

//
// SlidePartitionUp
//

template<typename T>
inline void
SlidePartitionUp
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2 )
{
#ifndef RELEASE
    CallStackEntry entry("SlidePartitionUp [Matrix]");
#endif
    View( AT, A0 );
    View2x1( AB, A1, A2 );
}

template<typename T,Distribution U,Distribution V>
inline void
SlidePartitionUp
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2 )
{
#ifndef RELEASE
    CallStackEntry entry("SlidePartitionUp [DistMatrix]");
#endif
    View( AT, A0 );
    View2x1( AB, A1, A2 );
}

template<typename T>
inline void
SlideLockedPartitionUp
( M& AT, const M& A0,
         const M& A1,
  M& AB, const M& A2 )
{
#ifndef RELEASE
    CallStackEntry entry("SlideLockedPartitionUp [Matrix]");
#endif
    LockedView( AT, A0 );
    LockedView2x1( AB, A1, A2 );
}

template<typename T,Distribution U,Distribution V>
inline void
SlideLockedPartitionUp
( DM& AT, const DM& A0,
          const DM& A1,
  DM& AB, const DM& A2 )
{
#ifndef RELEASE
    CallStackEntry entry("SlideLockedPartitionUp [DistMatrix]");
#endif
    LockedView( AT, A0 );
    LockedView2x1( AB, A1, A2 );
}

//
// SlidePartitionDown
//

template<typename T>
inline void
SlidePartitionDown
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2 )
{
#ifndef RELEASE
    CallStackEntry entry("SlidePartitionDown [Matrix]");
#endif
    View2x1( AT, A0, A1 );
    View( AB, A2 );
}

template<typename T,Distribution U,Distribution V>
inline void
SlidePartitionDown
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2 )
{
#ifndef RELEASE
    CallStackEntry entry("SlidePartitionDown [DistMatrix]");
#endif
    View2x1( AT, A0, A1 );
    View( AB, A2 );
}

template<typename T>
inline void
SlideLockedPartitionDown
( M& AT, const M& A0,
         const M& A1,
  M& AB, const M& A2 )
{
#ifndef RELEASE
    CallStackEntry entry("SlideLockedPartitionDown [Matrix]");
#endif
    LockedView2x1( AT, A0, A1 );
    LockedView( AB, A2 );
}

template<typename T,Distribution U,Distribution V>
inline void
SlideLockedPartitionDown
( DM& AT, const DM& A0,
          const DM& A1,
  DM& AB, const DM& A2 )
{
#ifndef RELEASE
    CallStackEntry entry("SlideLockedPartitionDown [DistMatrix]");
#endif
    LockedView2x1( AT, A0, A1 );
    LockedView( AB, A2 );
}

//
// SlidePartitionLeft
//

template<typename T>
inline void
SlidePartitionLeft
( M& AL, M& AR,
  M& A0, M& A1, M& A2 )
{
#ifndef RELEASE
    CallStackEntry entry("SlidePartitionLeft [Matrix]");
#endif
    View( AL, A0 );
    View1x2( AR, A1, A2 );
}

template<typename T,Distribution U,Distribution V>
inline void
SlidePartitionLeft
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2 )
{
#ifndef RELEASE
    CallStackEntry entry("SlidePartitionLeft [DistMatrix]");
#endif
    View( AL, A0 );
    View1x2( AR, A1, A2 );
}

template<typename T>
inline void
SlideLockedPartitionLeft
( M& AL, M& AR,
  const M& A0, const M& A1, const M& A2 )
{
#ifndef RELEASE
    CallStackEntry entry("SlideLockedPartitionLeft [Matrix]");
#endif
    LockedView( AL, A0 );
    LockedView1x2( AR, A1, A2 );
}

template<typename T,Distribution U,Distribution V>
inline void
SlideLockedPartitionLeft
( DM& AL, DM& AR,
  const DM& A0, const DM& A1, const DM& A2 )
{
#ifndef RELEASE
    CallStackEntry entry("SlideLockedPartitionLeft [DistMatrix]");
#endif
    LockedView( AL, A0 );
    LockedView1x2( AR, A1, A2 );
}

//
// SlidePartitionRight
//

template<typename T>
inline void
SlidePartitionRight
( M& AL, M& AR,
  M& A0, M& A1, M& A2 )
{
#ifndef RELEASE
    CallStackEntry entry("SlidePartitionRight [Matrix]");
#endif
    View1x2( AL, A0, A1 );
    View( AR, A2 );
}

template<typename T,Distribution U,Distribution V>
inline void
SlidePartitionRight
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2 )
{
#ifndef RELEASE
    CallStackEntry entry("SlidePartitionRight [DistMatrix]");
#endif
    View1x2( AL, A0, A1 );
    View( AR, A2 );
}

template<typename T>
inline void
SlideLockedPartitionRight
( M& AL, M& AR,
  const M& A0, const M& A1, const M& A2 )
{
#ifndef RELEASE
    CallStackEntry entry("SlideLockedPartitionRight [Matrix]");
#endif
    LockedView1x2( AL, A0, A1 );
    LockedView( AR, A2 );
}

template<typename T,Distribution U,Distribution V>
inline void
SlideLockedPartitionRight
( DM& AL, DM& AR,
  const DM& A0, const DM& A1, const DM& A2 )
{
#ifndef RELEASE
    CallStackEntry entry("SlideLockedPartitionRight [DistMatrix]");
#endif
    LockedView1x2( AL, A0, A1 );
    LockedView( AR, A2 );
}

//
// SlidePartitionUpDiagonal
//

template<typename T>
inline void
SlidePartitionUpDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22 )
{
#ifndef RELEASE
    CallStackEntry entry("SlidePartitionUpDiagonal [Matrix]");
#endif
    View( ATL, A00 );
    View1x2( ATR, A01, A02 );
    View2x1( ABL, A10, A20 );
    View2x2( ABR, A11, A12,
                  A21, A22 );
}

template<typename T,Distribution U,Distribution V>
inline void
SlidePartitionUpDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22 )
{
#ifndef RELEASE
    CallStackEntry entry("SlidePartitionUpDiagonal [DistMatrix]");
#endif
    View( ATL, A00 );
    View1x2( ATR, A01, A02 );
    View2x1( ABL, A10, A20 );
    View2x2( ABR, A11, A12,
                  A21, A22 );
}

template<typename T>
inline void
SlideLockedPartitionUpDiagonal
( M& ATL, M& ATR, const M& A00, const M& A01, const M& A02,
                  const M& A10, const M& A11, const M& A12,
  M& ABL, M& ABR, const M& A20, const M& A21, const M& A22 )
{
#ifndef RELEASE
    CallStackEntry entry("SlideLockedPartitionUpDiagonal [Matrix]");
#endif
    LockedView( ATL, A00 );
    LockedView1x2( ATR, A01, A02 );
    LockedView2x1( ABL, A10, A20 );
    LockedView2x2( ABR, A11, A12,
                        A21, A22 );
}

template<typename T,Distribution U,Distribution V>
inline void
SlideLockedPartitionUpDiagonal
( DM& ATL, DM& ATR, const DM& A00, const DM& A01, const DM& A02,
                    const DM& A10, const DM& A11, const DM& A12,
  DM& ABL, DM& ABR, const DM& A20, const DM& A21, const DM& A22 )
{
#ifndef RELEASE
    CallStackEntry entry("SlideLockedPartitionUpDiagonal [DistMatrix]");
#endif
    LockedView( ATL, A00 );
    LockedView1x2( ATR, A01, A02 );
    LockedView2x1( ABL, A10, A20 );
    LockedView2x2( ABR, A11, A12,
                        A21, A22 );
}

//
// SlidePartitionDownDiagonal
//

template<typename T>
inline void
SlidePartitionDownDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22 )
{
#ifndef RELEASE
    CallStackEntry entry("SlidePartitionDownDiagonal [Matrix]");
#endif
    View2x2( ATL, A00, A01,
                  A10, A11 );
    View2x1( ATR, A02, A12 );
    View1x2( ABL, A20, A21 );
    View( ABR, A22 );
}

template<typename T,Distribution U,Distribution V>
inline void
SlidePartitionDownDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22 )
{
#ifndef RELEASE
    CallStackEntry entry("SlidePartitionDownDiagonal [DistMatrix]");
#endif
    View2x2( ATL, A00, A01,
                  A10, A11 );
    View2x1( ATR, A02, A12 );
    View1x2( ABL, A20, A21 );
    View( ABR, A22 );
}

template<typename T>
inline void
SlideLockedPartitionDownDiagonal
( M& ATL, M& ATR, const M& A00, const M& A01, const M& A02,
                  const M& A10, const M& A11, const M& A12,
  M& ABL, M& ABR, const M& A20, const M& A21, const M& A22 )
{
#ifndef RELEASE
    CallStackEntry entry("SlideLockedPartitionDownDiagonal [Matrix]");
#endif
    LockedView2x2( ATL, A00, A01,
                        A10, A11 );
    LockedView2x1( ATR, A02, A12 );
    LockedView1x2( ABL, A20, A21 );
    LockedView( ABR, A22 );
}

template<typename T,Distribution U,Distribution V>
inline void
SlideLockedPartitionDownDiagonal
( DM& ATL, DM& ATR, const DM& A00, const DM& A01, const DM& A02,
                    const DM& A10, const DM& A11, const DM& A12,
  DM& ABL, DM& ABR, const DM& A20, const DM& A21, const DM& A22 )
{
#ifndef RELEASE
    CallStackEntry entry("SlideLockedPartitionDownDiagonal [DistMatrix]");
#endif
    LockedView2x2( ATL, A00, A01,
                        A10, A11 );
    LockedView2x1( ATR, A02, A12 );
    LockedView1x2( ABL, A20, A21 );
    LockedView( ABR, A22 );
}

#undef DM
#undef M

} // namespace elem

#endif // ifndef ELEM_CORE_SLIDEPARTITION_IMPL_HPP
