/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HILBERTSCHMIDT_HPP
#define LAPACK_HILBERTSCHMIDT_HPP

namespace elem {

// TODO: Think about using a more stable accumulation algorithm?

template<typename F> 
inline F
HilbertSchmidt( const Matrix<F>& A, const Matrix<F>& B )
{
#ifndef RELEASE
    CallStackEntry entry("HilbertSchmidt");
#endif
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        throw std::logic_error("Matrices must be the same size");
    F innerProd(0);
    const int width = A.Width();
    const int height = A.Height();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            innerProd += Conj(A.Get(i,j))*B.Get(i,j);
    return innerProd;
}

template<typename F,Distribution U,Distribution V> 
inline F
HilbertSchmidt( const DistMatrix<F,U,V>& A, const DistMatrix<F,U,V>& B )
{
#ifndef RELEASE
    CallStackEntry entry("HilbertSchmidt");
#endif
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        throw std::logic_error("Matrices must be the same size");
    if( A.Grid() != B.Grid() )
        throw std::logic_error("Grids must match");
    if( A.ColAlignment() != B.ColAlignment() || 
        A.RowAlignment() != B.RowAlignment() )
        throw std::logic_error("Matrices must be aligned");

    F localInnerProd(0);
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
            localInnerProd += Conj(A.GetLocal(iLoc,jLoc))*
                                   B.GetLocal(iLoc,jLoc);

    F innerProd;
    mpi::Comm comm = ReduceComm<U,V>( A.Grid() );
    mpi::AllReduce( &localInnerProd, &innerProd, 1, mpi::SUM, comm );
    return innerProd;
}

} // namespace elem

#endif // ifndef LAPACK_HILBERTSCHMIDT_HPP
