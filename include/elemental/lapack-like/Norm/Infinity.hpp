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

namespace elem {
namespace internal {

template<typename F> 
inline typename Base<F>::type
InfinityNorm( const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::InfinityNorm");
#endif
    typedef typename Base<F>::type R;

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

template<typename F,Distribution U,Distribution V> 
inline typename Base<F>::type
InfinityNorm( const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("internal::InfinityNorm");
#endif
    typedef typename Base<F>::type R;
    mpi::Comm colComm = ReduceColComm<U,V>( A.Grid() );
    mpi::Comm rowComm = ReduceRowComm<U,V>( A.Grid() );

    // Compute the partial row sums defined by our local matrix, A[U,V]
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
    mpi::AllReduce
    ( &myPartialRowSums[0], &myRowSums[0], localHeight, mpi::SUM, rowComm );

    // Find the maximum out of the row sums
    R myMaxRowSum = 0;
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
        myMaxRowSum = std::max( myMaxRowSum, myRowSums[iLocal] );

    // Find the global maximum row sum by searching over the U team
    R maxRowSum = 0;
    mpi::AllReduce( &myMaxRowSum, &maxRowSum, 1, mpi::MAX, colComm );
#ifndef RELEASE
    PopCallStack();
#endif
    return maxRowSum;
}

} // namespace internal
} // namespace elem

#endif // ifndef LAPACK_NORM_INFINITY_HPP
