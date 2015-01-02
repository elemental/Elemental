/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// TODO: Think about using a more stable accumulation algorithm?

template<typename T> 
T Dotu( const Matrix<T>& A, const Matrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Dotu"))
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Matrices must be the same size");
    T sum(0);
    const Int width = A.Width();
    const Int height = A.Height();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            sum += A.Get(i,j)*B.Get(i,j);
    return sum;
}

template<typename T> 
T Dotu( const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Dotu"))
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Matrices must be the same size");
    AssertSameGrids( A, B );
    if( A.DistData().colDist != B.DistData().colDist ||
        A.DistData().rowDist != B.DistData().rowDist )
        LogicError("Matrices must have the same distribution");
    if( A.ColAlign() != B.ColAlign() || 
        A.RowAlign() != B.RowAlign() )
        LogicError("Matrices must be aligned");

    T innerProd;
    if( A.Participating() )
    {
        T localInnerProd(0);
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                localInnerProd += A.GetLocal(iLoc,jLoc)*B.GetLocal(iLoc,jLoc);
        innerProd = mpi::AllReduce( localInnerProd, A.DistComm() );
    }
    mpi::Broadcast( innerProd, A.Root(), A.CrossComm() );
    return innerProd;
}

template<typename T>
T Dotu( const DistMultiVec<T>& A, const DistMultiVec<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Dotu"))
    if( !mpi::Congruent( A.Comm(), B.Comm() ) )
        LogicError("A and B must be congruent");
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("A and B must have the same dimensions");
    if( A.LocalHeight() != B.LocalHeight() )
        LogicError("A and B must have the same local heights");
    if( A.FirstLocalRow() != B.FirstLocalRow() )
        LogicError("A and B must own the same rows");

    T localInnerProd = 0;
    const Int localHeight = A.LocalHeight(); 
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )    
            localInnerProd += A.GetLocal(iLoc,j)*B.GetLocal(iLoc,j);
    return mpi::AllReduce( localInnerProd, A.Comm() );
}

#define PROTO(T) \
  template T Dotu( const Matrix<T>& A, const Matrix<T>& B ); \
  template T Dotu \
  ( const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B ); \
  template T Dotu( const DistMultiVec<T>& A, const DistMultiVec<T>& B );

#include "El/macros/Instantiate.h"

} // namespace El
