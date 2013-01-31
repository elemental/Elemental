/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_NORM_ENTRYWISE_HPP
#define LAPACK_NORM_ENTRYWISE_HPP

namespace elem {

template<typename F> 
inline typename Base<F>::type
EntrywiseNorm( const Matrix<F>& A, typename Base<F>::type p )
{
#ifndef RELEASE
    PushCallStack("EntrywiseNorm");
#endif
    // TODO: Make this more numerically stable
    typedef typename Base<F>::type R;
    R sum = 0;
    const int width = A.Width();
    const int height = A.Height();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            sum += Pow( Abs(A.Get(i,j)), p );
    const R norm = Pow( sum, 1/p );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F,Distribution U,Distribution V> 
inline typename Base<F>::type
EntrywiseNorm( const DistMatrix<F,U,V>& A, typename Base<F>::type p )
{
#ifndef RELEASE
    PushCallStack("EntrywiseNorm");
#endif
    typedef typename Base<F>::type R;
    R localSum = 0;
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            localSum += Pow( Abs(A.GetLocal(iLocal,jLocal)), p ); 

    R sum;
    mpi::Comm comm = ReduceComm<U,V>( A.Grid() );
    mpi::AllReduce( &localSum, &sum, 1, mpi::SUM, comm );
    const R norm = Pow( sum, 1/p );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem

#endif // ifndef LAPACK_NORM_ENTRYWISE_HPP
