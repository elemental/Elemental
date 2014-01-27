/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HILBERTSCHMIDT_HPP
#define ELEM_HILBERTSCHMIDT_HPP

namespace elem {

// TODO: Think about using a more stable accumulation algorithm?

template<typename F> 
inline F
HilbertSchmidt( const Matrix<F>& A, const Matrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("HilbertSchmidt"))
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Matrices must be the same size");
    F innerProd(0);
    const Int width = A.Width();
    const Int height = A.Height();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            innerProd += Conj(A.Get(i,j))*B.Get(i,j);
    return innerProd;
}

template<typename F,Dist U,Dist V> 
inline F
HilbertSchmidt( const DistMatrix<F,U,V>& A, const DistMatrix<F,U,V>& B )
{
    DEBUG_ONLY(CallStackEntry cse("HilbertSchmidt"))
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Matrices must be the same size");
    if( A.Grid() != B.Grid() )
        LogicError("Grids must match");
    if( A.ColAlign() != B.ColAlign() || A.RowAlign() != B.RowAlign() )
        LogicError("Matrices must be aligned");

    F innerProd;
    if( A.Participating() )
    {
        F localInnerProd(0);
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                localInnerProd += Conj(A.GetLocal(iLoc,jLoc))*
                                       B.GetLocal(iLoc,jLoc);
        innerProd = mpi::AllReduce( localInnerProd, A.DistComm() );
    }
    mpi::Broadcast( innerProd, A.Root(), A.CrossComm() );
    return innerProd;
}

} // namespace elem

#endif // ifndef ELEM_HILBERTSCHMIDT_HPP
