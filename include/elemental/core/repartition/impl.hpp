/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_REPARTITION_IMPL_HPP
#define ELEM_CORE_REPARTITION_IMPL_HPP

namespace elem {

// To make our life easier. Undef'd at the bottom of the header
#define M  Matrix<T>
#define DM DistMatrix<T,U,V>

//
// RepartitionUp
//

template<typename T>
inline void
RepartitionUp
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2, Int A1Height )
{
#ifndef RELEASE
    CallStackEntry cse("RepartitionUp [Matrix]");
#endif
    PartitionUp( AT, A0, A1, A1Height );
    View( A2, AB );
}

template<typename T,Distribution U,Distribution V>
inline void
RepartitionUp
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2, Int A1Height )
{
#ifndef RELEASE
    CallStackEntry cse("RepartitionUp [DistMatrix]");
#endif
    PartitionUp( AT, A0, A1, A1Height );
    View( A2, AB );
}

template<typename T>
inline void
LockedRepartitionUp
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2, Int A1Height )
{
#ifndef RELEASE
    CallStackEntry cse("LockedRepartitionUp [Matrix]");
#endif
    LockedPartitionUp( AT, A0, A1, A1Height );
    LockedView( A2, AB );
}

template<typename T,Distribution U,Distribution V>
inline void
LockedRepartitionUp
( const DM& AT, DM& A0,
                DM& A1,
  const DM& AB, DM& A2, Int A1Height )
{
#ifndef RELEASE
    CallStackEntry cse("LockedRepartitionUp [DistMatrix]");
#endif
    LockedPartitionUp( AT, A0, A1, A1Height );
    LockedView( A2, AB );
}

//
// RepartitionDown
//

template<typename T>
inline void
RepartitionDown
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2, Int A1Height )
{
#ifndef RELEASE
    CallStackEntry cse("RepartitionDown [Matrix]");
#endif
    View( A0, AT );
    PartitionDown( AB, A1, A2, A1Height );
}

template<typename T,Distribution U,Distribution V>
inline void
RepartitionDown
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2, Int A1Height )
{
#ifndef RELEASE
    CallStackEntry cse("RepartitionDown [DistMatrix]");
#endif
    View( A0, AT );
    PartitionDown( AB, A1, A2, A1Height );
}

template<typename T>
inline void
LockedRepartitionDown
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2, Int A1Height )
{
#ifndef RELEASE
    CallStackEntry cse("LockedRepartitionDown [Matrix]");
#endif
    LockedView( A0, AT );
    LockedPartitionDown( AB, A1, A2, A1Height );
}

template<typename T,Distribution U,Distribution V>
inline void
LockedRepartitionDown
( const DM& AT, DM& A0,
                DM& A1,
  const DM& AB, DM& A2, Int A1Height )
{
#ifndef RELEASE
    CallStackEntry cse("LockedRepartitionDown [DistMatrix]");
#endif
    LockedView( A0, AT );
    LockedPartitionDown( AB, A1, A2, A1Height );
}

//
// RepartitionLeft
//

template<typename T>
inline void
RepartitionLeft
( M& AL, M& AR,
  M& A0, M& A1, M& A2, Int A1Width )
{
#ifndef RELEASE
    CallStackEntry cse("RepartitionLeft [Matrix]");
#endif
    PartitionLeft( AL, A0, A1, A1Width );
    View( A2, AR );
}

template<typename T,Distribution U,Distribution V>
inline void
RepartitionLeft
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2, Int A1Width )
{
#ifndef RELEASE
    CallStackEntry cse("RepartitionLeft [DistMatrix]");
#endif
    PartitionLeft( AL, A0, A1, A1Width );
    View( A2, AR );
}

template<typename T>
inline void
LockedRepartitionLeft
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2, Int A1Width )
{
#ifndef RELEASE
    CallStackEntry cse("LockedRepartitionLeft [Matrix]");
#endif
    LockedPartitionLeft( AL, A0, A1, A1Width );
    LockedView( A2, AR );
}

template<typename T,Distribution U,Distribution V>
inline void
LockedRepartitionLeft
( const DM& AL, const DM& AR,
  DM& A0, DM& A1, DM& A2, Int A1Width )
{
#ifndef RELEASE
    CallStackEntry cse("LockedRepartitionLeft [DistMatrix]");
#endif
    LockedPartitionLeft( AL, A0, A1, A1Width );
    LockedView( A2, AR );
}

//
// RepartitionRight
//

template<typename T>
inline void
RepartitionRight
( M& AL, M& AR,
  M& A0, M& A1, M& A2, Int A1Width )
{
#ifndef RELEASE
    CallStackEntry cse("RepartitionRight [Matrix]");
#endif
    View( A0, AL );
    PartitionRight( AR, A1, A2, A1Width );
}

template<typename T,Distribution U,Distribution V>
inline void
RepartitionRight
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2, Int A1Width )
{
#ifndef RELEASE
    CallStackEntry cse("RepartitionRight [DistMatrix]");
#endif
    View( A0, AL );
    PartitionRight( AR, A1, A2, A1Width );
}

template<typename T>
inline void
LockedRepartitionRight
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2, Int A1Width )
{
#ifndef RELEASE
    CallStackEntry cse("LockedRepartitionRight [Matrix]");
#endif
    LockedView( A0, AL );
    LockedPartitionRight( AR, A1, A2, A1Width );
}

template<typename T,Distribution U,Distribution V>
inline void
LockedRepartitionRight
( const DM& AL, const DM& AR,
  DM& A0, DM& A1, DM& A2, Int A1Width )
{
#ifndef RELEASE
    CallStackEntry cse("LockedRepartitionRight [DistMatrix]");
#endif
    LockedView( A0, AL );
    LockedPartitionRight( AR, A1, A2, A1Width );
}

//
// RepartitionUpDiagonal
//

template<typename T>
inline void
RepartitionUpDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22, Int bsize )
{
#ifndef RELEASE
    CallStackEntry cse("RepartitionUpDiagonal [Matrix]");
#endif
    PartitionUpOffsetDiagonal
    ( ATL.Width()-ATL.Height(),
      ATL, A00, A01,
           A10, A11, bsize );
    PartitionUp( ATR, A02, A12, A11.Height() );
    PartitionLeft( ABL, A20, A21, A11.Width() );
    View( A22, ABR );
}

template<typename T,Distribution U,Distribution V>
inline void
RepartitionUpDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22, Int bsize )
{
#ifndef RELEASE
    CallStackEntry cse("RepartitionUpDiagonal [DistMatrix]");
#endif
    PartitionUpOffsetDiagonal
    ( ATL.Width()-ATL.Height(),
      ATL, A00, A01,
           A10, A11, bsize );
    PartitionUp( ATR, A02, A12, A11.Height() );
    PartitionLeft( ABL, A20, A21, A11.Width() );
    View( A22, ABR );
}

template<typename T>
inline void
LockedRepartitionUpDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22, Int bsize )
{
#ifndef RELEASE
    CallStackEntry cse("LockedRepartitionUpDiagonal [Matrix]");
#endif
    LockedPartitionUpOffsetDiagonal
    ( ATL.Width()-ATL.Height(),
      ATL, A00, A01,
           A10, A11, bsize );
    LockedPartitionUp( ATR, A02, A12, A11.Height() );
    LockedPartitionLeft( ABL, A20, A21, A11.Width() );
    LockedView( A22, ABR );
}

template<typename T,Distribution U,Distribution V>
inline void
LockedRepartitionUpDiagonal
( const DM& ATL, const DM& ATR, DM& A00, DM& A01, DM& A02,
                                DM& A10, DM& A11, DM& A12,
  const DM& ABL, const DM& ABR, DM& A20, DM& A21, DM& A22, Int bsize )
{
#ifndef RELEASE
    CallStackEntry cse("LockedRepartitionUpDiagonal [DistMatrix]");
#endif
    LockedPartitionUpOffsetDiagonal
    ( ATL.Width()-ATL.Height(),
      ATL, A00, A01,
           A10, A11, bsize );
    LockedPartitionUp( ATR, A02, A12, A11.Height() );
    LockedPartitionLeft( ABL, A20, A21, A11.Width() );
    LockedView( A22, ABR );
}

//
// RepartitionDownDiagonal
//

template<typename T>
inline void
RepartitionDownDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22, Int bsize )
{
#ifndef RELEASE
    CallStackEntry cse("RepartitionDownDiagonal [Matrix]");
#endif
    View( A00, ATL );
    PartitionDownDiagonal( ABR, A11, A12,
                                A21, A22, bsize );
    PartitionDown( ABL, A10, A20, A11.Height() );
    PartitionRight( ATR, A01, A02, A11.Width() );
}

template<typename T,Distribution U,Distribution V>
inline void
RepartitionDownDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22, Int bsize )
{
#ifndef RELEASE
    CallStackEntry cse("RepartitionDownDiagonal [DistMatrix]");
#endif
    View( A00, ATL );
    PartitionDownDiagonal( ABR, A11, A12,
                                A21, A22, bsize );
    PartitionDown( ABL, A10, A20, A11.Height() );
    PartitionRight( ATR, A01, A02, A11.Width() );
}

template<typename T>
inline void
LockedRepartitionDownDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22, Int bsize )
{
#ifndef RELEASE
    CallStackEntry cse("LockedRepartitionDownDiagonal [Matrix]");
#endif
    LockedView( A00, ATL );
    LockedPartitionDownDiagonal( ABR, A11, A12,
                                      A21, A22, bsize );
    LockedPartitionDown( ABL, A10, A20, A11.Height() );
    LockedPartitionRight( ATR, A01, A02, A11.Width() );
}

template<typename T,Distribution U,Distribution V>
inline void
LockedRepartitionDownDiagonal
( const DM& ATL, const DM& ATR, DM& A00, DM& A01, DM& A02,
                                DM& A10, DM& A11, DM& A12,
  const DM& ABL, const DM& ABR, DM& A20, DM& A21, DM& A22, Int bsize )
{
#ifndef RELEASE
    CallStackEntry cse("LockedRepartitionDownDiagonal [DistMatrix]");
#endif
    LockedView( A00, ATL );
    LockedPartitionDownDiagonal( ABR, A11, A12,
                                      A21, A22, bsize );
    LockedPartitionDown( ABL, A10, A20, A11.Height() );
    LockedPartitionRight( ATR, A01, A02, A11.Width() );
}

#undef DM
#undef M

} // namespace elem

#endif // ifndef ELEM_CORE_REPARTITION_IMPL_HPP
