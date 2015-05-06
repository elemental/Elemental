/*
   Copyright (c) 2009-2015, Jack Poulson
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
    typedef Base<F> R;
    R maxRowSum = 0;
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int i=0; i<height; ++i )
    {
        R rowSum = 0;
        for( Int j=0; j<width; ++j )
            rowSum += Abs(A.Get(i,j));
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
        vector<Real> myPartialRowSums( localHeight );
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            myPartialRowSums[iLoc] = 0;
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                myPartialRowSums[iLoc] += Abs(A.GetLocal(iLoc,jLoc));
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
  ( UpperOrLower uplo, const AbstractDistMatrix<T>& A ); 

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
