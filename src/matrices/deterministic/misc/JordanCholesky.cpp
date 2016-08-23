/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>
#include <El/matrices.hpp>

namespace El {

template<typename T> 
void JordanCholesky( Matrix<T>& A, Int n )
{
    DEBUG_CSE
    Zeros( A, n, n );

    // Set the main diagonal equal to five everywhere but the top-left entry
    FillDiagonal( A, T(5) );
    if( n > 0 )
        A.Set( 0, 0, T(1) );

    // Set the rest of the tridiagonal to 2
    FillDiagonal( A, T(2), 1 );
    FillDiagonal( A, T(2), -1 );
}

template<typename T>
void JordanCholesky( AbstractDistMatrix<T>& A, Int n )
{
    DEBUG_CSE
    Zeros( A, n, n );

    // Set the main diagonal equal to five everywhere but the top-left entry
    FillDiagonal( A, T(5) );
    if( n > 0 )
        A.Set( 0, 0, T(1) );

    // Set the rest of the tridiagonal to 2
    FillDiagonal( A, T(2), 1 );
    FillDiagonal( A, T(2), -1 );
}

template<typename T>
void JordanCholesky( SparseMatrix<T>& A, Int n )
{
    DEBUG_CSE
    Zeros( A, n, n );
    A.Reserve( 3*n );
    
    for( Int e=0; e<n; ++e )
    {
        if( e == 0 )
            A.QueueUpdate( e, e, T(1) );
        else
            A.QueueUpdate( e, e, T(5) );
        if( e > 0 )
            A.QueueUpdate( e, e-1, T(2) );
        if( e < n-1 )
            A.QueueUpdate( e, e+1, T(2) );
    }
    A.ProcessQueues();
}

template<typename T>
void JordanCholesky( DistSparseMatrix<T>& A, Int n )
{
    DEBUG_CSE
    Zeros( A, n, n );

    const Int localHeight = A.LocalHeight();
    A.Reserve( 3*localHeight );

    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = A.GlobalRow(iLoc);
        if( i == 0 )
            A.QueueUpdate( i, i, T(1) );
        else
            A.QueueUpdate( i, i, T(5) );

        if( i > 0 )
            A.QueueUpdate( i, i-1, T(2) );
        if( i < n-1 )
            A.QueueUpdate( i, i+1, T(2) );
    }
    A.ProcessQueues();
}

#define PROTO(T) \
  template void JordanCholesky( Matrix<T>& A, Int n ); \
  template void JordanCholesky( AbstractDistMatrix<T>& A, Int n ); \
  template void JordanCholesky( SparseMatrix<T>& A, Int n ); \
  template void JordanCholesky( DistSparseMatrix<T>& A, Int n );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
