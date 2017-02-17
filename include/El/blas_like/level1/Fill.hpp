/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_FILL_HPP
#define EL_BLAS_FILL_HPP

namespace El {

template<typename T>
void Fill( Matrix<T>& A, T alpha )
{
    EL_DEBUG_CSE
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            A(i,j) = alpha;
}

template<typename T>
void Fill( AbstractDistMatrix<T>& A, T alpha )
{
    EL_DEBUG_CSE
    Fill( A.Matrix(), alpha );
}

template<typename T>
void Fill( DistMultiVec<T>& A, T alpha )
{
    EL_DEBUG_CSE
    Fill( A.Matrix(), alpha );
}

template<typename T>
void Fill( SparseMatrix<T>& A, T alpha )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    A.Resize( m, n );
    Zero( A );
    if( alpha != T(0) )
    {
        A.Reserve( m*n );
        for( Int i=0; i<m; ++i )
            for( Int j=0; j<n; ++j )
                A.QueueUpdate( i, j, alpha );
        A.ProcessQueues();
    }
}

template<typename T>
void Fill( DistSparseMatrix<T>& A, T alpha )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    A.Resize( m, n );
    Zero( A );
    if( alpha != T(0) )
    {
        const Int localHeight = A.LocalHeight();
        A.Reserve( localHeight*n );
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            for( Int j=0; j<n; ++j )
                A.QueueLocalUpdate( iLoc, j, alpha );
        A.ProcessLocalQueues();
    }
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void Fill( Matrix<T>& A, T alpha ); \
  EL_EXTERN template void Fill( AbstractDistMatrix<T>& A, T alpha ); \
  EL_EXTERN template void Fill( DistMultiVec<T>& A, T alpha ); \
  EL_EXTERN template void Fill( SparseMatrix<T>& A, T alpha ); \
  EL_EXTERN template void Fill( DistSparseMatrix<T>& A, T alpha );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_FILL_HPP
