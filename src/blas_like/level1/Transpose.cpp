/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void Transpose( const Matrix<T>& A, Matrix<T>& B, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Transpose"))
    const Int m = A.Height();
    const Int n = A.Width();
    B.Resize( n, m );
    // TODO: Optimize this routine
    const T* ABuf = A.LockedBuffer();
          T* BBuf = B.Buffer();
    const Int ldA = A.LDim();
    const Int ldB = B.LDim();
    if( conjugate )
    {
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<m; ++i )
                BBuf[j+i*ldB] = Conj(ABuf[i+j*ldA]);
    }
    else
    {
        copy::util::InterleaveMatrix
        ( m, n, 
          ABuf, 1,   ldA,
          BBuf, ldB, 1 );
    }
}

template<typename T>
void Transpose
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Transpose"))
    const DistData ADistData = A.DistData();
    const DistData BDistData = B.DistData();
    if( ADistData.colDist == BDistData.rowDist &&
        ADistData.rowDist == BDistData.colDist &&
        ((ADistData.colAlign==BDistData.rowAlign) || !B.RowConstrained()) &&
        ((ADistData.rowAlign==BDistData.colAlign) || !B.ColConstrained()) )
    {
        B.Align( A.RowAlign(), A.ColAlign() );
        B.Resize( A.Width(), A.Height() );
        Transpose( A.LockedMatrix(), B.Matrix(), conjugate );
    }
    else
    {
        std::unique_ptr<AbstractDistMatrix<T>> 
            C( B.ConstructTranspose(A.Grid(),A.Root()) );
        C->AlignRowsWith( B.DistData() );
        C->AlignColsWith( B.DistData() );
        Copy( A, *C );
        B.Resize( A.Width(), A.Height() );
        Transpose( C->LockedMatrix(), B.Matrix(), conjugate );
    }
}

template<typename T>
void Transpose
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B, 
  bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Transpose"))
    const BlockDistData ADistData = A.DistData();
    const BlockDistData BDistData = B.DistData();
    if( ADistData.colDist == BDistData.rowDist &&
        ADistData.rowDist == BDistData.colDist &&
        ((ADistData.colAlign    == BDistData.rowAlign && 
          ADistData.blockHeight == BDistData.blockWidth &&
          ADistData.colCut      == BDistData.rowCut) || !B.RowConstrained()) &&
        ((ADistData.rowAlign   == BDistData.colAlign && 
          ADistData.blockWidth == BDistData.blockHeight &&
          ADistData.rowCut     == BDistData.colCut) || !B.ColConstrained()))
    {
        B.Align
        ( A.BlockWidth(), A.BlockHeight(), 
          A.RowAlign(), A.ColAlign(), A.RowCut(), A.ColCut() );
        B.Resize( A.Width(), A.Height() );
        Transpose( A.LockedMatrix(), B.Matrix(), conjugate );
    }
    else
    {
        std::unique_ptr<AbstractBlockDistMatrix<T>> 
            C( B.ConstructTranspose(A.Grid(),A.Root()) );
        C->AlignRowsWith( B.DistData() );
        C->AlignColsWith( B.DistData() );
        Copy( A, *C );
        B.Resize( A.Width(), A.Height() );
        Transpose( C->LockedMatrix(), B.Matrix(), conjugate );
    }
}

template<typename T>
void Transpose
( const SparseMatrix<T>& A, SparseMatrix<T>& B, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Transpose"))
    Zeros( B, A.Width(), A.Height() );
    TransposeAxpy( T(1), A, B, conjugate );
}

template<typename T>
void Transpose
( const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Transpose"))
    B.SetComm( A.Comm() );
    Zeros( B, A.Width(), A.Height() );
    TransposeAxpy( T(1), A, B, conjugate );
}

#define PROTO(T) \
  template void Transpose( const Matrix<T>& A, Matrix<T>& B, bool conjugate ); \
  template void Transpose \
  ( const AbstractDistMatrix<T>& A, \
          AbstractDistMatrix<T>& B, bool conjugate ); \
  template void Transpose \
  ( const AbstractBlockDistMatrix<T>& A, \
          AbstractBlockDistMatrix<T>& B, bool conjugate ); \
  template void Transpose \
  ( const SparseMatrix<T>& A, \
          SparseMatrix<T>& B, bool conjugate ); \
  template void Transpose \
  ( const DistSparseMatrix<T>& A, \
          DistSparseMatrix<T>& B, bool conjugate );

#include "El/macros/Instantiate.h"

} // namespace El
