/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_NORM_INFINITY_HPP
#define LAPACK_NORM_INFINITY_HPP

#include "elemental/lapack-like/Norm/One.hpp"

namespace elem {

template<typename F> 
inline BASE(F)
InfinityNorm( const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("InfinityNorm");
#endif
    typedef BASE(F) R;
    R maxRowSum = 0;
    const int height = A.Height();
    const int width = A.Width();
    for( int i=0; i<height; ++i )
    {
        R rowSum = 0;
        for( int j=0; j<width; ++j )
            rowSum += Abs(A.Get(i,j));
        maxRowSum = std::max( maxRowSum, rowSum );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return maxRowSum;
}

template<typename F>
inline BASE(F)
HermitianInfinityNorm( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("HermitianInfinityNorm");
#endif
    typedef BASE(F) R;
    R maxRowSum = HermitianOneNorm( uplo, A );
#ifndef RELEASE
    PopCallStack();
#endif
    return maxRowSum;
}

template<typename F>
inline BASE(F)
SymmetricInfinityNorm( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("SymmetricInfinityNorm");
#endif
    typedef BASE(F) R;
    const R norm = HermitianInfinityNorm( uplo, A );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F,Distribution U,Distribution V> 
inline BASE(F)
InfinityNorm( const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("InfinityNorm");
#endif
    // Compute the partial row sums defined by our local matrix, A[U,V]
    typedef BASE(F) R;
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    std::vector<R> myPartialRowSums( localHeight );
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
    {
        myPartialRowSums[iLocal] = 0;
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
            myPartialRowSums[iLocal] += Abs(A.GetLocal(iLocal,jLocal));
    }

    // Sum our partial row sums to get the row sums over A[U,* ]
    std::vector<R> myRowSums( localHeight );
    mpi::Comm rowComm = ReduceRowComm<U,V>( A.Grid() );
    mpi::AllReduce
    ( &myPartialRowSums[0], &myRowSums[0], localHeight, mpi::SUM, rowComm );

    // Find the maximum out of the row sums
    R myMaxRowSum = 0;
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
        myMaxRowSum = std::max( myMaxRowSum, myRowSums[iLocal] );

    // Find the global maximum row sum by searching over the U team
    R maxRowSum = 0;
    mpi::Comm colComm = ReduceColComm<U,V>( A.Grid() );
    mpi::AllReduce( &myMaxRowSum, &maxRowSum, 1, mpi::MAX, colComm );
#ifndef RELEASE
    PopCallStack();
#endif
    return maxRowSum;
}

template<typename F>
inline BASE(F)
HermitianInfinityNorm
( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("HermitianInfinityNorm");
#endif
    typedef BASE(F) R;
    R maxRowSum = HermitianOneNorm( uplo, A );
#ifndef RELEASE
    PopCallStack();
#endif
    return maxRowSum;
}

template<typename F>
inline BASE(F)
SymmetricInfinityNorm( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("SymmetricInfinityNorm");
#endif
    typedef BASE(F) R;
    const R norm = HermitianInfinityNorm( uplo, A );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem

#endif // ifndef LAPACK_NORM_INFINITY_HPP
