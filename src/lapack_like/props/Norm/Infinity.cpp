/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F> 
Base<F> InfinityNorm( const Matrix<F>& A )
{
    DEBUG_ONLY(CSE cse("InfinityNorm"))
    typedef Base<F> Real;
    const Int height = A.Height();
    const Int width = A.Width();
    const F* ABuf = A.LockedBuffer();
    const Int ALDim = A.LDim();

    Real maxRowSum = 0;
    for( Int i=0; i<height; ++i )
    {
        Real rowSum = 0;
        for( Int j=0; j<width; ++j )
            rowSum += Abs(ABuf[i+j*ALDim]);
        maxRowSum = Max( maxRowSum, rowSum );
    }
    return maxRowSum;
}

template<typename F>
Base<F> HermitianInfinityNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CSE cse("HermitianInfinityNorm"))
    return HermitianOneNorm( uplo, A );
}

template<typename F>
Base<F> SymmetricInfinityNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CSE cse("SymmetricInfinityNorm"))
    return HermitianInfinityNorm( uplo, A );
}

template<typename F> 
Base<F> InfinityNorm( const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("InfinityNorm"))
    // Compute the partial row sums defined by our local matrix, A[U,V]
    typedef Base<F> Real;


    Real norm;
    if( A.Participating() )
    {
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        const F* ABuf = A.LockedBuffer();
        const Int ALDim = A.LDim();

        vector<Real> myPartialRowSums( localHeight );
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            myPartialRowSums[iLoc] = 0;
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                myPartialRowSums[iLoc] += Abs(ABuf[iLoc+jLoc*ALDim]);
        }

        // Sum our partial row sums to get the row sums over A[U,* ]
        vector<Real> myRowSums( localHeight );
        mpi::AllReduce
        ( myPartialRowSums.data(), myRowSums.data(), localHeight, A.RowComm() );

        // Find the maximum out of the row sums
        Real myMaxRowSum = 0;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            myMaxRowSum = Max( myMaxRowSum, myRowSums[iLoc] );

        // Find the global maximum row sum by searching over the U team
        norm = mpi::AllReduce( myMaxRowSum, mpi::MAX, A.ColComm() );
    }
    mpi::Broadcast( norm, A.Root(), A.CrossComm() );
    return norm;
}

template<typename F>
Base<F> HermitianInfinityNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("HermitianInfinityNorm"))
    return HermitianOneNorm( uplo, A );
}

template<typename F>
Base<F> SymmetricInfinityNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("SymmetricInfinityNorm"))
    return HermitianInfinityNorm( uplo, A );
}

template<typename F> 
Base<F> InfinityNorm( const SparseMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("InfinityNorm"))
    typedef Base<F> Real;
    const Int height = A.Height();
    const F* valBuf = A.LockedValueBuffer();
    const Int* offsetBuf = A.LockedOffsetBuffer();

    Real maxRowSum = 0;
    for( Int i=0; i<height; ++i ) 
    {
        const Int thisOff = offsetBuf[i];
        const Int nextOff = offsetBuf[i+1];
        Real rowSum = 0;
        for( Int k=thisOff; k<nextOff; ++k )
            rowSum += Abs(valBuf[k]);
        maxRowSum = Max(rowSum,maxRowSum);
    }

    return maxRowSum;
}

template<typename F> 
Base<F> InfinityNorm( const DistSparseMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("InfinityNorm"))
    typedef Base<F> Real;
    const Int localHeight = A.LocalHeight();
    const F* valBuf = A.LockedValueBuffer();
    const Int* offsetBuf = A.LockedOffsetBuffer();

    Real maxLocRowSum = 0;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc ) 
    {
        const Int thisOff = offsetBuf[iLoc];
        const Int nextOff = offsetBuf[iLoc+1];
        Real rowSum = 0;
        for( Int k=thisOff; k<nextOff; ++k )
            rowSum += Abs(valBuf[k]);
        maxLocRowSum = Max(rowSum,maxLocRowSum);
    }

    return mpi::AllReduce( maxLocRowSum, mpi::MAX, A.Comm() );
}

#define PROTO(T) \
  template Base<T> InfinityNorm( const Matrix<T>& A ); \
  template Base<T> InfinityNorm ( const AbstractDistMatrix<T>& A ); \
  template Base<T> HermitianInfinityNorm \
  ( UpperOrLower uplo, const Matrix<T>& A ); \
  template Base<T> HermitianInfinityNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<T>& A ); \
  template Base<T> SymmetricInfinityNorm \
  ( UpperOrLower uplo, const Matrix<T>& A ); \
  template Base<T> SymmetricInfinityNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<T>& A ); \
  template Base<T> InfinityNorm( const SparseMatrix<T>& A ); \
  template Base<T> InfinityNorm( const DistSparseMatrix<T>& A );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
