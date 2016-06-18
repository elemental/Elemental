/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_RESHAPE_HPP
#define EL_BLAS_RESHAPE_HPP

namespace El {

template<typename T>
void Reshape
(       Int mNew,
        Int nNew, 
  const Matrix<T>& A,
        Matrix<T>& B )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    if( m*n != mNew*nNew )
        LogicError
        ("Reshape from ",m," x ",n," to ",mNew," x ",nNew,
         " did not preserve the total number of entries");

    B.Resize( mNew, nNew );
    Zero( B );
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            const Int iNew = (i+j*m) % mNew;
            const Int jNew = (i+j*m) / mNew;
            B.Set( iNew, jNew, A.Get(i,j) );
        }
    }
}

template<typename T>
Matrix<T> Reshape( Int mNew, Int nNew, const Matrix<T>& A )
{
    DEBUG_CSE
    Matrix<T> B;
    Reshape( mNew, nNew, A, B );
    return B;
}

// TODO: Merge with implementation of GetSubmatrix via a function which maps
//       the coordinates in A to the coordinates in B
template<typename T>
void Reshape
(       Int mNew,
        Int nNew, 
  const AbstractDistMatrix<T>& A,
        AbstractDistMatrix<T>& B )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    const Grid& g = A.Grid();
    if( m*n != mNew*nNew )
        LogicError
        ("Reshape from ",m," x ",n," to ",mNew," x ",nNew,
         " did not preserve the total number of entries");

    B.SetGrid( g ); 
    B.Resize( mNew, nNew );
    Zero( B );
    
    B.Reserve( mLocal*nLocal );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<mLocal; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            const Int iNew = (i+j*m) % mNew;
            const Int jNew = (i+j*m) / mNew;
            B.QueueUpdate( iNew, jNew, A.GetLocal(iLoc,jLoc) );
        }
    }
    B.ProcessQueues();
}

template<typename T>
DistMatrix<T> Reshape( Int mNew, Int nNew, const AbstractDistMatrix<T>& A )
{
    DistMatrix<T> B( A.Grid() );
    Reshape( mNew, nNew, A, B );
    return B;
}

template<typename T>
void Reshape
(       Int mNew,
        Int nNew,
  const SparseMatrix<T>& A,
        SparseMatrix<T>& B )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numEntries = A.NumEntries();
    if( m*n != mNew*nNew )
        LogicError
        ("Reshape from ",m," x ",n," to ",mNew," x ",nNew,
         " did not preserve the total number of entries");

    B.Resize( mNew, nNew );
    Zero( B );
    B.Reserve( numEntries );

    // Insert the nonzeros
    for( Int e=0; e<numEntries; ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        const Int iNew = (i+j*m) % mNew;
        const Int jNew = (i+j*m) / mNew;
        B.QueueUpdate( iNew, jNew, A.Value(e) );
    }
    B.ProcessQueues();
}

template<typename T>
SparseMatrix<T> Reshape( Int mNew, Int nNew, const SparseMatrix<T>& A )
{
    SparseMatrix<T> B;
    Reshape( mNew, nNew, A, B );
    return B;
}

template<typename T>
void Reshape
(       Int mNew,
        Int nNew, 
  const DistSparseMatrix<T>& A,
        DistSparseMatrix<T>& B )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numEntries = A.NumLocalEntries();
    if( m*n != mNew*nNew )
        LogicError
        ("Reshape from ",m," x ",n," to ",mNew," x ",nNew,
         " did not preserve the total number of entries");

    B.SetComm( A.Comm() );
    B.Resize( mNew, nNew );
    Zero( B );

    B.Reserve( numEntries );
    for( Int e=0; e<numEntries; ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        const Int iNew = (i+j*m) % mNew;
        const Int jNew = (i+j*m) / mNew;
        B.QueueUpdate( iNew, jNew, A.Value(e) );
    }
    B.ProcessQueues();
}

template<typename T>
DistSparseMatrix<T> Reshape( Int mNew, Int nNew, const DistSparseMatrix<T>& A )
{
    DistSparseMatrix<T> B(A.Comm());
    Reshape( mNew, nNew, A, B );
    return B;
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void Reshape \
  (       Int mNew, \
          Int nNew, \
    const Matrix<T>& A, \
          Matrix<T>& B ); \
  EL_EXTERN template Matrix<T> Reshape \
  ( Int mNew, Int nNew, const Matrix<T>& A ); \
  EL_EXTERN template void Reshape \
  (       Int mNew, \
          Int nNew, \
    const AbstractDistMatrix<T>& A, \
          AbstractDistMatrix<T>& B ); \
  EL_EXTERN template DistMatrix<T> Reshape \
  ( Int mNew, Int nNew, const AbstractDistMatrix<T>& A ); \
  EL_EXTERN template void Reshape \
  (       Int mNew, \
          Int nNew, \
    const SparseMatrix<T>& A, \
          SparseMatrix<T>& B ); \
  EL_EXTERN template SparseMatrix<T> Reshape \
  ( Int mNew, Int nNew, const SparseMatrix<T>& A ); \
  EL_EXTERN template void Reshape \
  (       Int mNew, \
          Int nNew, \
    const DistSparseMatrix<T>& A, \
          DistSparseMatrix<T>& B ); \
  EL_EXTERN template DistSparseMatrix<T> Reshape \
  ( Int mNew, Int nNew, const DistSparseMatrix<T>& A );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_RESHAPE_HPP
