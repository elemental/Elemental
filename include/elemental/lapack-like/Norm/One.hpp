/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_NORM_ONE_HPP
#define LAPACK_NORM_ONE_HPP

namespace elem {
namespace internal {

template<typename F>
inline typename Base<F>::type
OneNorm( const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::OneNorm");
#endif
    typedef typename Base<F>::type R;

    R maxColSum = 0;
    const int height = A.Height();
    const int width = A.Width();
    for( int j=0; j<width; ++j )
    {
        R colSum = 0;
        for( int i=0; i<height; ++i )
            colSum += Abs(A.Get(i,j));
        maxColSum = std::max( maxColSum, colSum );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return maxColSum;
}

template<typename F,Distribution U,Distribution V>
inline typename Base<F>::type
OneNorm( const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("internal::OneNorm");
#endif
    typedef typename Base<F>::type R;

    mpi::Comm colComm = ReduceColComm<U,V>( A.Grid() );
    mpi::Comm rowComm = ReduceRowComm<U,V>( A.Grid() );

    // Compute the partial column sums defined by our local matrix, A[U,V]
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    std::vector<R> myPartialColSums( localWidth );
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        myPartialColSums[jLocal] = 0;
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            myPartialColSums[jLocal] += Abs(A.GetLocal(iLocal,jLocal));
    }

    // Sum our partial column sums to get the column sums over A[* ,V]
    std::vector<R> myColSums( localWidth );
    mpi::AllReduce
    ( &myPartialColSums[0], &myColSums[0], localWidth, mpi::SUM, colComm );

    // Find the maximum out of the column sums
    R myMaxColSum = 0;
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
        myMaxColSum = std::max( myMaxColSum, myColSums[jLocal] );

    // Find the global maximum column sum by searching over the MR team
    R maxColSum = 0;
    mpi::AllReduce( &myMaxColSum, &maxColSum, 1, mpi::MAX, rowComm );
#ifndef RELEASE
    PopCallStack();
#endif
    return maxColSum;
}

} // namespace internal
} // namespace elem

#endif // ifndef LAPACK_NORM_ONE_HPP
