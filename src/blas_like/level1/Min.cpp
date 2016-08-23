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

template<typename Real,typename>
Real Min( const Matrix<Real>& A )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Real* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();

    Real value = limits::Max<Real>();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            value = Min(value,ABuf[i+j*ALDim]);
    return value;
}

template<typename Real,typename>
Real Min( const AbstractDistMatrix<Real>& A )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( !A.Grid().InGrid() )
          LogicError("Viewing processes are not allowed");
    )
    Real value = limits::Max<Real>();
    if( A.Participating() )
    {
        // Store the index/value of the local pivot candidate
        const Int mLocal = A.LocalHeight();
        const Int nLocal = A.LocalWidth();
        const Real* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                value = Min(value,ABuf[iLoc+jLoc*ALDim]);

        value = mpi::AllReduce( value, mpi::MIN, A.DistComm() );
    }
    mpi::Broadcast( value, A.Root(), A.CrossComm() );
    return value;
}

template<typename Real,typename>
Real SymmetricMin( UpperOrLower uplo, const Matrix<Real>& A )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
    )
    const Int n = A.Width();
    const Real* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();

    Real value = limits::Max<Real>();
    if( uplo == LOWER )
    {
        for( Int j=0; j<n; ++j )
            for( Int i=j; i<n; ++i )
                value = Min(value,ABuf[i+j*ALDim]);
    }
    else
    {
        for( Int j=0; j<n; ++j ) 
            for( Int i=0; i<=j; ++i )
                value = Min(value,ABuf[i+j*ALDim]);
    }
    return value;
}

template<typename Real,typename>
Real SymmetricMin( UpperOrLower uplo, const AbstractDistMatrix<Real>& A )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( !A.Grid().InGrid() )
          LogicError("Viewing processes are not allowed");
    )

    Real value = limits::Max<Real>();
    if( A.Participating() )
    {
        const Int mLocal = A.LocalHeight();
        const Int nLocal = A.LocalWidth();
        const Real* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();
        if( uplo == LOWER )
        {
            for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int mLocBefore = A.LocalRowOffset(j);
                for( Int iLoc=mLocBefore; iLoc<mLocal; ++iLoc )
                    value = Min(value,ABuf[iLoc+jLoc*ALDim]);
            }
        }
        else
        {
            for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                const Int mLocBefore = A.LocalRowOffset(j+1);
                for( Int iLoc=0; iLoc<mLocBefore; ++iLoc )
                    value = Min(value,ABuf[iLoc+jLoc*ALDim]);
            }
        }
        value = mpi::AllReduce( value, mpi::MIN, A.DistComm() );
    }
    mpi::Broadcast( value, A.Root(), A.CrossComm() );
    return value;
}

#define PROTO(Real) \
  template Real Min( const Matrix<Real>& x ); \
  template Real Min( const AbstractDistMatrix<Real>& x ); \
  template Real SymmetricMin( UpperOrLower uplo, const Matrix<Real>& A ); \
  template Real SymmetricMin \
  ( UpperOrLower uplo, const AbstractDistMatrix<Real>& A );

#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
