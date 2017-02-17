/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>

namespace El {

template<typename Ring>
Base<Ring> MaxAbs( const Matrix<Ring>& A )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Ring* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();

    Base<Ring> value = 0;
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            value = Max(value,Abs(ABuf[i+j*ALDim]));
    return value;
}

template<typename Ring>
Base<Ring> MaxAbs( const AbstractDistMatrix<Ring>& A )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( !A.Grid().InGrid() )
          LogicError("Viewing processes are not allowed");
    )
    Base<Ring> value = 0;
    if( A.Participating() )
    {
        // Store the index/value of the local pivot candidate
        const Int mLocal = A.LocalHeight();
        const Int nLocal = A.LocalWidth();
        const Ring* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                value = Max(value,Abs(ABuf[iLoc+jLoc*ALDim]));

        value = mpi::AllReduce( value, mpi::MAX, A.DistComm() );
    }
    mpi::Broadcast( value, A.Root(), A.CrossComm() );
    return value;
}

template<typename Ring>
Base<Ring> SymmetricMaxAbs( UpperOrLower uplo, const Matrix<Ring>& A )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
    )
    const Int n = A.Width();
    const Ring* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();

    Base<Ring> value = 0;
    if( uplo == LOWER )
    {
        for( Int j=0; j<n; ++j )
            for( Int i=j; i<n; ++i )
                value = Max(value,Abs(ABuf[i+j*ALDim]));
    }
    else
    {
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<=j; ++i )
                value = Max(value,Abs(ABuf[i+j*ALDim]));
    }
    return value;
}

template<typename Ring>
Base<Ring> SymmetricMaxAbs
( UpperOrLower uplo, const AbstractDistMatrix<Ring>& A )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( !A.Grid().InGrid() )
          LogicError("Viewing processes are not allowed");
    )

    Base<Ring> value = 0;
    if( A.Participating() )
    {
        const Int mLocal = A.LocalHeight();
        const Int nLocal = A.LocalWidth();
        const Ring* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();
        if( uplo == LOWER )
        {
            for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int mLocBefore = A.LocalRowOffset(j);
                for( Int iLoc=mLocBefore; iLoc<mLocal; ++iLoc )
                    value = Max(value,Abs(ABuf[iLoc+jLoc*ALDim]));
            }
        }
        else
        {
            for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int mLocBefore = A.LocalRowOffset(j+1);
                for( Int iLoc=0; iLoc<mLocBefore; ++iLoc )
                    value = Max(value,Abs(ABuf[iLoc+jLoc*ALDim]));
            }
        }
        value = mpi::AllReduce( value, mpi::MAX, A.DistComm() );
    }
    mpi::Broadcast( value, A.Root(), A.CrossComm() );
    return value;
}

#define PROTO(Ring) \
  template Base<Ring> MaxAbs( const Matrix<Ring>& x ); \
  template Base<Ring> MaxAbs( const AbstractDistMatrix<Ring>& x ); \
  template Base<Ring> SymmetricMaxAbs \
  ( UpperOrLower uplo, const Matrix<Ring>& x ); \
  template Base<Ring> SymmetricMaxAbs \
  ( UpperOrLower uplo, const AbstractDistMatrix<Ring>& x );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
