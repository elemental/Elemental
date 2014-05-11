/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_NORM_INFINITY_HPP
#define ELEM_NORM_INFINITY_HPP

#include ELEM_ONENORM_INC

namespace elem {

template<typename F> 
inline Base<F>
InfinityNorm( const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("InfinityNorm"))
    typedef Base<F> R;
    R maxRowSum = 0;
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int i=0; i<height; ++i )
    {
        R rowSum = 0;
        for( Int j=0; j<width; ++j )
            rowSum += Abs(A.Get(i,j));
        maxRowSum = std::max( maxRowSum, rowSum );
    }
    return maxRowSum;
}

template<typename F>
inline Base<F>
HermitianInfinityNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianInfinityNorm"))
    return HermitianOneNorm( uplo, A );
}

template<typename F>
inline Base<F>
SymmetricInfinityNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricInfinityNorm"))
    return HermitianInfinityNorm( uplo, A );
}

template<typename F,Dist U,Dist V> 
inline Base<F>
InfinityNorm( const DistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("InfinityNorm"))
    // Compute the partial row sums defined by our local matrix, A[U,V]
    typedef Base<F> Real;
    Real norm;
    if( A.Participating() )
    {
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        std::vector<Real> myPartialRowSums( localHeight );
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            myPartialRowSums[iLoc] = 0;
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                myPartialRowSums[iLoc] += Abs(A.GetLocal(iLoc,jLoc));
        }

        // Sum our partial row sums to get the row sums over A[U,* ]
        std::vector<Real> myRowSums( localHeight );
        mpi::AllReduce
        ( myPartialRowSums.data(), myRowSums.data(), localHeight, A.RowComm() );

        // Find the maximum out of the row sums
        Real myMaxRowSum = 0;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            myMaxRowSum = std::max( myMaxRowSum, myRowSums[iLoc] );

        // Find the global maximum row sum by searching over the U team
        norm = mpi::AllReduce( myMaxRowSum, mpi::MAX, A.ColComm() );
    }
    mpi::Broadcast( norm, A.Root(), A.CrossComm() );
    return norm;
}

template<typename F,Dist U,Dist V>
inline Base<F>
HermitianInfinityNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianInfinityNorm"))
    return HermitianOneNorm( uplo, A );
}

template<typename F,Dist U,Dist V>
inline Base<F>
SymmetricInfinityNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricInfinityNorm"))
    return HermitianInfinityNorm( uplo, A );
}

} // namespace elem

#endif // ifndef ELEM_NORM_INFINITY_HPP
