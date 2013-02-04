/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_NORM_ENTRYWISEONE_HPP
#define LAPACK_NORM_ENTRYWISEONE_HPP

namespace elem {

template<typename F> 
inline typename Base<F>::type
EntrywiseOneNorm( const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("EntrywiseOneNorm");
#endif
    typedef typename Base<F>::type R;
    R norm = 0;
    const int width = A.Width();
    const int height = A.Height();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            norm += Abs(A.Get(i,j));
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F,Distribution U,Distribution V> 
inline typename Base<F>::type
EntrywiseOneNorm( const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("EntrywiseOneNorm");
#endif
    typedef typename Base<F>::type R;
    R localSum = 0;
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            localSum += Abs(A.GetLocal(iLocal,jLocal)); 

    R norm;
    mpi::Comm comm = ReduceComm<U,V>( A.Grid() );
    mpi::AllReduce( &localSum, &norm, 1, mpi::SUM, comm );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem

#endif // ifndef LAPACK_NORM_ENTRYWISEONE_HPP
