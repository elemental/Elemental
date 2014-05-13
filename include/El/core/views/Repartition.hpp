/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_VIEWS_REPARTITION_HPP
#define EL_VIEWS_REPARTITION_HPP

#include "./Partition.hpp"

namespace El {

// To make our life easier. Undef'd at the bottom of the header
#define M  Matrix<T>
#define DM DistMatrix<T,U,V>

// Repartition upwards from the bottom
// ===================================

template<typename T>
inline void
RepartitionUp
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2, Int A1Height=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("RepartitionUp"))
    PartitionUp( AT, A0, A1, A1Height );
    View( A2, AB );
}

template<typename T,Dist U,Dist V>
inline void
RepartitionUp
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2, Int A1Height=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("RepartitionUp"))
    PartitionUp( AT, A0, A1, A1Height );
    View( A2, AB );
}

template<typename T>
inline void
LockedRepartitionUp
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2, Int A1Height=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedRepartitionUp"))
    LockedPartitionUp( AT, A0, A1, A1Height );
    LockedView( A2, AB );
}

template<typename T,Dist U,Dist V>
inline void
LockedRepartitionUp
( const DM& AT, DM& A0,
                DM& A1,
  const DM& AB, DM& A2, Int A1Height=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedRepartitionUp"))
    LockedPartitionUp( AT, A0, A1, A1Height );
    LockedView( A2, AB );
}

// Repartition downwards from the top
// ==================================

template<typename T>
inline void
RepartitionDown
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2, Int A1Height=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("RepartitionDown"))
    View( A0, AT );
    PartitionDown( AB, A1, A2, A1Height );
}

template<typename T,Dist U,Dist V>
inline void
RepartitionDown
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2, Int A1Height=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("RepartitionDown"))
    View( A0, AT );
    PartitionDown( AB, A1, A2, A1Height );
}

template<typename T>
inline void
LockedRepartitionDown
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2, Int A1Height=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedRepartitionDown"))
    LockedView( A0, AT );
    LockedPartitionDown( AB, A1, A2, A1Height );
}

template<typename T,Dist U,Dist V>
inline void
LockedRepartitionDown
( const DM& AT, DM& A0,
                DM& A1,
  const DM& AB, DM& A2, Int A1Height=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedRepartitionDown"))
    LockedView( A0, AT );
    LockedPartitionDown( AB, A1, A2, A1Height );
}

// Repartition leftwards from the right
// ====================================

template<typename T>
inline void
RepartitionLeft
( M& AL, M& AR,
  M& A0, M& A1, M& A2, Int A1Width=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("RepartitionLeft"))
    PartitionLeft( AL, A0, A1, A1Width );
    View( A2, AR );
}

template<typename T,Dist U,Dist V>
inline void
RepartitionLeft
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2, Int A1Width=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("RepartitionLeft"))
    PartitionLeft( AL, A0, A1, A1Width );
    View( A2, AR );
}

template<typename T>
inline void
LockedRepartitionLeft
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2, Int A1Width=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedRepartitionLeft"))
    LockedPartitionLeft( AL, A0, A1, A1Width );
    LockedView( A2, AR );
}

template<typename T,Dist U,Dist V>
inline void
LockedRepartitionLeft
( const DM& AL, const DM& AR,
  DM& A0, DM& A1, DM& A2, Int A1Width=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedRepartitionLeft"))
    LockedPartitionLeft( AL, A0, A1, A1Width );
    LockedView( A2, AR );
}

// Repartition rightward from the left
// ===================================

template<typename T>
inline void
RepartitionRight
( M& AL, M& AR,
  M& A0, M& A1, M& A2, Int A1Width=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("RepartitionRight"))
    View( A0, AL );
    PartitionRight( AR, A1, A2, A1Width );
}

template<typename T,Dist U,Dist V>
inline void
RepartitionRight
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2, Int A1Width=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("RepartitionRight"))
    View( A0, AL );
    PartitionRight( AR, A1, A2, A1Width );
}

template<typename T>
inline void
LockedRepartitionRight
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2, Int A1Width=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedRepartitionRight"))
    LockedView( A0, AL );
    LockedPartitionRight( AR, A1, A2, A1Width );
}

template<typename T,Dist U,Dist V>
inline void
LockedRepartitionRight
( const DM& AL, const DM& AR,
  DM& A0, DM& A1, DM& A2, Int A1Width=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedRepartitionRight"))
    LockedView( A0, AL );
    LockedPartitionRight( AR, A1, A2, A1Width );
}

// Repartition upwards on a diagonal
// =================================

template<typename T>
inline void
RepartitionUpDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("RepartitionUpDiagonal"))
    PartitionUpOffsetDiagonal
    ( ATL.Width()-ATL.Height(),
      ATL, A00, A01,
           A10, A11, bsize );
    PartitionUp( ATR, A02, A12, A11.Height() );
    PartitionLeft( ABL, A20, A21, A11.Width() );
    View( A22, ABR );
}

template<typename T,Dist U,Dist V>
inline void
RepartitionUpDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22, Int bsize=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("RepartitionUpDiagonal"))
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
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedRepartitionUpDiagonal"))
    LockedPartitionUpOffsetDiagonal
    ( ATL.Width()-ATL.Height(),
      ATL, A00, A01,
           A10, A11, bsize );
    LockedPartitionUp( ATR, A02, A12, A11.Height() );
    LockedPartitionLeft( ABL, A20, A21, A11.Width() );
    LockedView( A22, ABR );
}

template<typename T,Dist U,Dist V>
inline void
LockedRepartitionUpDiagonal
( const DM& ATL, const DM& ATR, DM& A00, DM& A01, DM& A02,
                                DM& A10, DM& A11, DM& A12,
  const DM& ABL, const DM& ABR, DM& A20, DM& A21, DM& A22, 
  Int bsize=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedRepartitionUpDiagonal"))
    LockedPartitionUpOffsetDiagonal
    ( ATL.Width()-ATL.Height(),
      ATL, A00, A01,
           A10, A11, bsize );
    LockedPartitionUp( ATR, A02, A12, A11.Height() );
    LockedPartitionLeft( ABL, A20, A21, A11.Width() );
    LockedView( A22, ABR );
}

// Repartition downwards on a diagonal
// ===================================

template<typename T>
inline void
RepartitionDownDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("RepartitionDownDiagonal"))
    View( A00, ATL );
    PartitionDownDiagonal( ABR, A11, A12,
                                A21, A22, bsize );
    PartitionDown( ABL, A10, A20, A11.Height() );
    PartitionRight( ATR, A01, A02, A11.Width() );
}

template<typename T,Dist U,Dist V>
inline void
RepartitionDownDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22, Int bsize=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("RepartitionDownDiagonal"))
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
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedRepartitionDownDiagonal"))
    LockedView( A00, ATL );
    LockedPartitionDownDiagonal( ABR, A11, A12,
                                      A21, A22, bsize );
    LockedPartitionDown( ABL, A10, A20, A11.Height() );
    LockedPartitionRight( ATR, A01, A02, A11.Width() );
}

template<typename T,Dist U,Dist V>
inline void
LockedRepartitionDownDiagonal
( const DM& ATL, const DM& ATR, DM& A00, DM& A01, DM& A02,
                                DM& A10, DM& A11, DM& A12,
  const DM& ABL, const DM& ABR, DM& A20, DM& A21, DM& A22, 
  Int bsize=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedRepartitionDownDiagonal"))
    LockedView( A00, ATL );
    LockedPartitionDownDiagonal( ABR, A11, A12,
                                      A21, A22, bsize );
    LockedPartitionDown( ABL, A10, A20, A11.Height() );
    LockedPartitionRight( ATR, A01, A02, A11.Width() );
}

#undef DM
#undef M

} // namespace El

#endif // ifndef EL_VIEWS_REPARTITION_HPP
