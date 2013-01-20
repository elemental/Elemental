/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_DOTU_HPP
#define LAPACK_DOTU_HPP

namespace elem {

// TODO: Think about using a more stable accumulation algorithm?

template<typename F> 
inline F
Dotu( const Matrix<F>& A, const Matrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("Dotu");
#endif
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        throw std::logic_error("Matrices must be the same size");
    F sum(0);
    const int width = A.Width();
    const int height = A.Height();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            sum += A.Get(i,j)*B.Get(i,j);
#ifndef RELEASE
    PopCallStack();
#endif
    return sum;
}

template<typename F,Distribution U,Distribution V> 
inline F
Dotu( const DistMatrix<F,U,V>& A, const DistMatrix<F,U,V>& B )
{
#ifndef RELEASE
    PushCallStack("Dotu");
#endif
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        throw std::logic_error("Matrices must be the same size");
    if( A.Grid() != B.Grid() )
        throw std::logic_error("Grids must match");
    if( A.ColAlignment() != B.ColAlignment() || 
        A.RowAlignment() != B.RowAlignment() )
        throw std::logic_error("Matrices must be aligned");

    F localSum(0);
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            localSum += A.GetLocal(iLocal,jLocal)*
                        B.GetLocal(iLocal,jLocal);

    F sum;
    mpi::Comm comm = ReduceComm<U,V>( A.Grid() );
    mpi::AllReduce( &localSum, &sum, 1, mpi::SUM, comm );
#ifndef RELEASE
    PopCallStack();
#endif
    return sum;
}

} // namespace elem

#endif // ifndef LAPACK_HILBERTSCHMIDT_HPP
