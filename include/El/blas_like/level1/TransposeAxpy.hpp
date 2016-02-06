/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_TRANSPOSEAXPY_HPP
#define EL_BLAS_TRANSPOSEAXPY_HPP

namespace El {

template<typename T,typename S>
void TransposeAxpy
(       S alphaS,
  const Matrix<T>& X,
        Matrix<T>& Y,
        bool conjugate )
{
    DEBUG_ONLY(CSE cse("TransposeAxpy"))
    const T alpha = T(alphaS);
    const Int mX = X.Height();
    const Int nX = X.Width();
    const Int nY = Y.Width();
    const Int ldX = X.LDim();
    const Int ldY = Y.LDim();
    const T* XBuf = X.LockedBuffer();
          T* YBuf = Y.Buffer();
    // If X and Y are vectors, we can allow one to be a column and the other
    // to be a row. Otherwise we force X and Y to be the same dimension.
    if( mX == 1 || nX == 1 )
    {
        const Int lengthX = ( nX==1 ? mX : nX );
        const Int incX = ( nX==1 ? 1  : ldX );
        const Int incY = ( nY==1 ? 1  : ldY );
        DEBUG_ONLY(
          const Int mY = Y.Height();
          const Int lengthY = ( nY==1 ? mY : nY );
          if( lengthX != lengthY )
              LogicError("Nonconformal TransposeAxpy");
        )
        if( conjugate )
            for( Int j=0; j<lengthX; ++j )
                YBuf[j*incY] += alpha*Conj(XBuf[j*incX]);
        else
            blas::Axpy( lengthX, alpha, XBuf, incX, YBuf, incY );
    }
    else
    {
        DEBUG_ONLY(
          const Int mY = Y.Height();
          if( mX != nY || nX != mY )
              LogicError("Nonconformal TransposeAxpy");
        )
        if( nX <= mX )
        {
            if( conjugate )
                for( Int j=0; j<nX; ++j )
                    for( Int i=0; i<mX; ++i )
                        YBuf[j+i*ldY] += alpha*Conj(XBuf[i+j*ldX]);
            else
                for( Int j=0; j<nX; ++j )
                    blas::Axpy( mX, alpha, &XBuf[j*ldX], 1, &YBuf[j], ldY );
        }
        else
        {
            if( conjugate )
                for( Int i=0; i<mX; ++i )
                    for( Int j=0; j<nX; ++j )
                        YBuf[j+i*ldY] += alpha*Conj(XBuf[i+j*ldX]);
            else
                for( Int i=0; i<mX; ++i )
                    blas::Axpy( nX, alpha, &XBuf[i], ldX, &YBuf[i*ldY], 1 );
        }
    }
}

template<typename T,typename S>
void TransposeAxpy
(       S alphaS,
  const SparseMatrix<T>& X,
        SparseMatrix<T>& Y,
        bool conjugate )
{
    DEBUG_ONLY(CSE cse("TransposeAxpy"))
    if( X.Height() != Y.Width() || X.Width() != Y.Height() )
        LogicError("X and Y must have transposed dimensions");
    const T alpha = T(alphaS);
    const Int numEntries = X.NumEntries();
    Y.Reserve( Y.NumEntries()+numEntries );
    for( Int k=0; k<numEntries; ++k ) 
    {
        const T value = alpha*( conjugate ? Conj(X.Value(k)) : X.Value(k) );
        Y.QueueUpdate( X.Col(k), X.Row(k), value );
    }
    Y.ProcessQueues();
}

template<typename T,typename S>
void TransposeAxpy
(       S alphaS,
  const ElementalMatrix<T>& A,
        ElementalMatrix<T>& B,
        bool conjugate )
{
    DEBUG_ONLY(
      CSE cse("TransposeAxpy");
      AssertSameGrids( A, B );
      if( A.Height() != B.Width() || A.Width() != B.Height() )
          LogicError("A and B must have transposed dimensions");
    )
    const T alpha = T(alphaS);

    const ElementalData ADistData = A.DistData();
    const ElementalData BDistData = B.DistData();
    if( ADistData.colDist == BDistData.rowDist &&
        ADistData.rowDist == BDistData.colDist &&
        ADistData.colAlign==BDistData.rowAlign &&
        ADistData.rowAlign==BDistData.colAlign )
    {
        TransposeAxpy( alpha, A.LockedMatrix(), B.Matrix(), conjugate );
    }
    else
    {
        unique_ptr<ElementalMatrix<T>>
            C( B.ConstructTranspose(A.Grid(),A.Root()) );
        C->AlignRowsWith( B.DistData() );
        C->AlignColsWith( B.DistData() );
        Copy( A, *C );
        TransposeAxpy( alpha, C->LockedMatrix(), B.Matrix(), conjugate );
    }
}

template<typename T,typename S>
void TransposeAxpy
(       S alphaS,
  const DistSparseMatrix<T>& A,
        DistSparseMatrix<T>& B,
        bool conjugate )
{
    DEBUG_ONLY(CSE cse("TransposeAxpy"))
    if( A.Height() != B.Width() || A.Width() != B.Height() )
        LogicError("A and B must have transposed dimensions");
    if( A.Comm() != B.Comm() )
        LogicError("A and B must have the same communicator");
    
    const Int numLocalEntries = A.NumLocalEntries();

    T alpha(alphaS);
    B.Reserve( B.NumLocalEntries()+numLocalEntries, numLocalEntries );
    for( Int e=0; e<numLocalEntries; ++e )
        B.QueueUpdate
        ( A.Col(e), A.Row(e),
          alpha*(conjugate ? Conj(A.Value(e)) : A.Value(e)) );
    B.ProcessQueues();
}

template<typename T,typename S>
void AdjointAxpy( S alphaS, const Matrix<T>& X, Matrix<T>& Y )
{
    DEBUG_ONLY(CSE cse("AdjointAxpy"))
    TransposeAxpy( alphaS, X, Y, true );
}

template<typename T,typename S>
void AdjointAxpy( S alphaS, const SparseMatrix<T>& X, SparseMatrix<T>& Y )
{
    DEBUG_ONLY(CSE cse("AdjointAxpy"))
    TransposeAxpy( alphaS, X, Y, true );
}

template<typename T,typename S>
void AdjointAxpy
( S alphaS, const ElementalMatrix<T>& X, ElementalMatrix<T>& Y )
{
    DEBUG_ONLY(CSE cse("AdjointAxpy"))
    TransposeAxpy( alphaS, X, Y, true );
}

template<typename T,typename S>
void AdjointAxpy
( S alphaS, const DistSparseMatrix<T>& X, DistSparseMatrix<T>& Y )
{
    DEBUG_ONLY(CSE cse("AdjointAxpy"))
    TransposeAxpy( alphaS, X, Y, true );
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO_TYPES(T,S) \
  EL_EXTERN template void TransposeAxpy \
  (       S alpha, \
    const Matrix<T>& A, \
          Matrix<T>& B, \
          bool conjugate ); \
  EL_EXTERN template void TransposeAxpy \
  (       S alpha, \
    const ElementalMatrix<T>& A, \
          ElementalMatrix<T>& B, \
          bool conjugate ); \
  EL_EXTERN template void TransposeAxpy \
  (       S alpha, \
    const SparseMatrix<T>& A, \
          SparseMatrix<T>& B, \
          bool conjugate ); \
  EL_EXTERN template void TransposeAxpy \
  (       S alpha, \
    const DistSparseMatrix<T>& A, \
          DistSparseMatrix<T>& B, \
          bool conjugate ); \
  EL_EXTERN template void AdjointAxpy \
  (       S alpha, \
    const Matrix<T>& A, \
          Matrix<T>& B ); \
  EL_EXTERN template void AdjointAxpy \
  (       S alpha, \
    const ElementalMatrix<T>& A, \
          ElementalMatrix<T>& B ); \
  EL_EXTERN template void AdjointAxpy \
  (       S alpha, \
    const SparseMatrix<T>& A, \
          SparseMatrix<T>& B ); \
  EL_EXTERN template void AdjointAxpy \
  (       S alpha, \
    const DistSparseMatrix<T>& A, \
          DistSparseMatrix<T>& B );

#define PROTO_INT(T) PROTO_TYPES(T,T)

#define PROTO_REAL(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,T)

#define PROTO_COMPLEX(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,Base<T>) \
  PROTO_TYPES(T,T)

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

#undef PROTO_TYPES
#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_TRANSPOSEAXPY_HPP
