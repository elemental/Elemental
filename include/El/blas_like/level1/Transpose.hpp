/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_TRANSPOSE_HPP
#define EL_BLAS_TRANSPOSE_HPP

namespace El {

namespace transpose {

template<typename T>
void ColFilter
( const ElementalMatrix<T>& A,
        ElementalMatrix<T>& B, bool conjugate );
template<typename T>
void ColFilter
( const BlockMatrix<T>& A,
        BlockMatrix<T>& B, bool conjugate );

template<typename T>
void RowFilter
( const ElementalMatrix<T>& A,
        ElementalMatrix<T>& B, bool conjugate );
template<typename T>
void RowFilter
( const BlockMatrix<T>& A,
        BlockMatrix<T>& B, bool conjugate );

template<typename T>
void PartialColFilter
( const ElementalMatrix<T>& A,
        ElementalMatrix<T>& B, bool conjugate );
template<typename T>
void PartialColFilter
( const BlockMatrix<T>& A,
        BlockMatrix<T>& B, bool conjugate );

template<typename T>
void PartialRowFilter
( const ElementalMatrix<T>& A,
        ElementalMatrix<T>& B, bool conjugate );
template<typename T>
void PartialRowFilter
( const BlockMatrix<T>& A,
        BlockMatrix<T>& B, bool conjugate );

template<typename T>
void ColAllGather
( const ElementalMatrix<T>& A,
        ElementalMatrix<T>& B, bool conjugate );
template<typename T>
void ColAllGather
( const BlockMatrix<T>& A,
        BlockMatrix<T>& B, bool conjugate );

template<typename T>
void PartialColAllGather
( const ElementalMatrix<T>& A,
        ElementalMatrix<T>& B, bool conjugate );
template<typename T>
void PartialColAllGather
( const BlockMatrix<T>& A,
        BlockMatrix<T>& B, bool conjugate );

} // namespace transpose

template<typename T>
void Transpose( const Matrix<T>& A, Matrix<T>& B, bool conjugate )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    B.Resize( n, m );
#ifdef EL_HAVE_MKL
    Orientation orient = ( conjugate ? ADJOINT : TRANSPOSE );
    mkl::omatcopy
    ( orient, m, n, T(1), A.LockedBuffer(), A.LDim(), B.Buffer(), B.LDim() );
#else
    // OpenBLAS's {i,o}matcopy routines where disabled for the reasons detailed
    // in src/core/imports/openblas.cpp

    // TODO(poulson): Optimize this routine
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
#endif
}

template<typename T>
void Transpose
( const ElementalMatrix<T>& A,
        ElementalMatrix<T>& B,
  bool conjugate )
{
    EL_DEBUG_CSE
    const auto AData = A.DistData();
    const auto BData = B.DistData();

    // NOTE: The following are ordered in terms of increasing cost
    if( AData.colDist == BData.rowDist &&
        AData.rowDist == BData.colDist &&
        ((AData.colAlign==BData.rowAlign) || !B.RowConstrained()) &&
        ((AData.rowAlign==BData.colAlign) || !B.ColConstrained()) )
    {
        B.Align( A.RowAlign(), A.ColAlign() );
        B.Resize( A.Width(), A.Height() );
        Transpose( A.LockedMatrix(), B.Matrix(), conjugate );
    }
    else if( AData.colDist == BData.rowDist &&
             AData.rowDist == Collect(BData.colDist) )
    {
        transpose::ColFilter( A, B, conjugate );
    }
    else if( AData.colDist == Collect(BData.rowDist) &&
             AData.rowDist == BData.colDist )
    {
        transpose::RowFilter( A, B, conjugate );
    }
    else if( AData.colDist == BData.rowDist &&
             AData.rowDist == Partial(BData.colDist) )
    {
        transpose::PartialColFilter( A, B, conjugate );
    }
    else if( AData.colDist == Partial(BData.rowDist) &&
             AData.rowDist == BData.colDist )
    {
        transpose::PartialRowFilter( A, B, conjugate );
    }
    else if( Partial(AData.colDist) == BData.rowDist &&
             AData.rowDist          == BData.colDist )
    {
        transpose::PartialColAllGather( A, B, conjugate );
    }
    else if( Collect(AData.colDist) == BData.rowDist &&
             AData.rowDist          == BData.colDist )
    {
        transpose::ColAllGather( A, B, conjugate );
    }
    else
    {
        unique_ptr<ElementalMatrix<T>>
            C( B.ConstructTranspose(A.Grid(),A.Root()) );
        C->AlignWith( BData );
        Copy( A, *C );
        B.Resize( A.Width(), A.Height() );
        Transpose( C->LockedMatrix(), B.Matrix(), conjugate );
    }
}

template<typename T>
void Transpose
( const BlockMatrix<T>& A,
        BlockMatrix<T>& B,
  bool conjugate )
{
    EL_DEBUG_CSE
    const auto AData = A.DistData();
    const auto BData = B.DistData();
    if( AData.colDist == BData.rowDist &&
        AData.rowDist == BData.colDist &&
        ((AData.colAlign    == BData.rowAlign &&
          AData.blockHeight == BData.blockWidth &&
          AData.colCut      == BData.rowCut) || !B.RowConstrained()) &&
        ((AData.rowAlign   == BData.colAlign &&
          AData.blockWidth == BData.blockHeight &&
          AData.rowCut     == BData.colCut) || !B.ColConstrained()))
    {
        B.Align
        ( A.BlockWidth(), A.BlockHeight(),
          A.RowAlign(), A.ColAlign(), A.RowCut(), A.ColCut() );
        B.Resize( A.Width(), A.Height() );
        Transpose( A.LockedMatrix(), B.Matrix(), conjugate );
    }
    else if( AData.colDist == BData.rowDist &&
             AData.rowDist == Collect(BData.colDist) )
    {
        transpose::ColFilter( A, B, conjugate );
    }
    else if( AData.colDist == Collect(BData.rowDist) &&
             AData.rowDist == BData.colDist )
    {
        transpose::RowFilter( A, B, conjugate );
    }
    else if( AData.colDist == BData.rowDist &&
             AData.rowDist == Partial(BData.colDist) )
    {
        transpose::PartialColFilter( A, B, conjugate );
    }
    else if( AData.colDist == Partial(BData.rowDist) &&
             AData.rowDist == BData.colDist )
    {
        transpose::PartialRowFilter( A, B, conjugate );
    }
    else if( Partial(AData.colDist) == BData.rowDist &&
             AData.rowDist          == BData.colDist )
    {
        transpose::PartialColAllGather( A, B, conjugate );
    }
    else if( Collect(AData.colDist) == BData.rowDist &&
             AData.rowDist          == BData.colDist )
    {
        transpose::ColAllGather( A, B, conjugate );
    }
    else
    {
        unique_ptr<BlockMatrix<T>>
            C( B.ConstructTranspose(A.Grid(),A.Root()) );
        C->AlignWith( BData );
        Copy( A, *C );
        B.Resize( A.Width(), A.Height() );
        Transpose( C->LockedMatrix(), B.Matrix(), conjugate );
    }
}

template<typename T>
void Transpose
( const AbstractDistMatrix<T>& A,
        AbstractDistMatrix<T>& B,
  bool conjugate )
{
    EL_DEBUG_CSE
    if( A.Wrap() == ELEMENT && B.Wrap() == ELEMENT )
    {
        const auto& ACast = static_cast<const ElementalMatrix<T>&>(A);
              auto& BCast = static_cast<      ElementalMatrix<T>&>(B);
        Transpose( ACast, BCast, conjugate );
    }
    else if( A.Wrap() == BLOCK  && B.Wrap() == BLOCK )
    {
        const auto& ACast = static_cast<const BlockMatrix<T>&>(A);
              auto& BCast = static_cast<      BlockMatrix<T>&>(B);
        Transpose( ACast, BCast, conjugate );
    }
    else if( A.Wrap() == ELEMENT ) // && B.Wrap() == BLOCK
    {
        auto& BCast = static_cast<BlockMatrix<T>&>(B);
        unique_ptr<BlockMatrix<T>>
            C( BCast.ConstructTranspose(A.Grid(),A.Root()) );
        C->AlignWith( BCast );
        Copy( A, *C );
        BCast.Resize( A.Width(), A.Height() );
        Transpose( C->LockedMatrix(), BCast.Matrix(), conjugate );
    }
    else  // A.Wrap() == BLOCK && B.Wrap() == ELEMENT
    {
        auto& BCast = static_cast<ElementalMatrix<T>&>(B);
        unique_ptr<ElementalMatrix<T>>
            C( BCast.ConstructTranspose(A.Grid(),A.Root()) );
        C->AlignWith( BCast );
        Copy( A, *C );
        BCast.Resize( A.Width(), A.Height() );
        Transpose( C->LockedMatrix(), BCast.Matrix(), conjugate );
    }
}

template<typename T>
void Transpose
( const SparseMatrix<T>& A, SparseMatrix<T>& B, bool conjugate )
{
    EL_DEBUG_CSE
    B.Resize( A.Width(), A.Height() );
    Zero( B, false );
    TransposeAxpy( T(1), A, B, conjugate );
}

template<typename T>
void Transpose
( const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B, bool conjugate )
{
    EL_DEBUG_CSE
    B.SetGrid( A.Grid() );
    B.Resize( A.Width(), A.Height() );
    Zero( B, false );
    TransposeAxpy( T(1), A, B, conjugate );
}

template<typename T>
void Adjoint( const Matrix<T>& A, Matrix<T>& B )
{
    EL_DEBUG_CSE
    Transpose( A, B, true );
}

template<typename T>
void Adjoint( const ElementalMatrix<T>& A, ElementalMatrix<T>& B )
{
    EL_DEBUG_CSE
    Transpose( A, B, true );
}

template<typename T>
void Adjoint
( const BlockMatrix<T>& A, BlockMatrix<T>& B )
{
    EL_DEBUG_CSE
    Transpose( A, B, true );
}

template<typename T>
void Adjoint
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B )
{
    EL_DEBUG_CSE
    Transpose( A, B, true );
}

template<typename T>
void Adjoint( const SparseMatrix<T>& A, SparseMatrix<T>& B )
{
    EL_DEBUG_CSE
    Transpose( A, B, true );
}

template<typename T>
void Adjoint( const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B )
{
    EL_DEBUG_CSE
    Transpose( A, B, true );
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void Transpose \
  ( const Matrix<T>& A, Matrix<T>& B, bool conjugate ); \
  EL_EXTERN template void Transpose \
  ( const ElementalMatrix<T>& A, ElementalMatrix<T>& B, bool conjugate ); \
  EL_EXTERN template void Transpose \
  ( const BlockMatrix<T>& A, BlockMatrix<T>& B, bool conjugate ); \
  EL_EXTERN template void Transpose \
  ( const AbstractDistMatrix<T>& A, \
          AbstractDistMatrix<T>& B, bool conjugate ); \
  EL_EXTERN template void Transpose \
  ( const SparseMatrix<T>& A, SparseMatrix<T>& B, bool conjugate ); \
  EL_EXTERN template void Transpose \
  ( const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B, bool conjugate ); \
  EL_EXTERN template void Adjoint \
  ( const Matrix<T>& A, Matrix<T>& B ); \
  EL_EXTERN template void Adjoint \
  ( const ElementalMatrix<T>& A, ElementalMatrix<T>& B ); \
  EL_EXTERN template void Adjoint \
  ( const BlockMatrix<T>& A, BlockMatrix<T>& B ); \
  EL_EXTERN template void Adjoint \
  ( const AbstractDistMatrix<T>& A, \
          AbstractDistMatrix<T>& B ); \
  EL_EXTERN template void Adjoint \
  ( const SparseMatrix<T>& A, SparseMatrix<T>& B ); \
  EL_EXTERN template void Adjoint \
  ( const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#include <El/blas_like/level1/Transpose/ColAllGather.hpp>
#include <El/blas_like/level1/Transpose/ColFilter.hpp>
#include <El/blas_like/level1/Transpose/PartialColAllGather.hpp>
#include <El/blas_like/level1/Transpose/PartialColFilter.hpp>
#include <El/blas_like/level1/Transpose/PartialRowFilter.hpp>
#include <El/blas_like/level1/Transpose/RowFilter.hpp>

#endif // ifndef EL_BLAS_TRANSPOSE_HPP
