/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_TRANSPOSECONTRACT_HPP
#define EL_BLAS_TRANSPOSECONTRACT_HPP

namespace El {

template<typename T>
void TransposeContract
( const ElementalMatrix<T>& A,
        ElementalMatrix<T>& B, bool conjugate )
{
    DEBUG_CSE
    const Dist U = B.ColDist();
    const Dist V = B.RowDist();
    if( A.ColDist() == V && A.RowDist() == Partial(U) )
    {
        Transpose( A, B, conjugate );
    }
    else
    {
        unique_ptr<ElementalMatrix<T>> 
            ASumFilt( B.ConstructTranspose(B.Grid(),B.Root()) );
        if( B.ColConstrained() )
            ASumFilt->AlignRowsWith( B, true );
        if( B.RowConstrained() )
            ASumFilt->AlignColsWith( B, true );
        Contract( A, *ASumFilt );
        if( !B.ColConstrained() )
            B.AlignColsWith( *ASumFilt, false );
        if( !B.RowConstrained() )
            B.AlignRowsWith( *ASumFilt, false );
        // We should have ensured that the alignments match
        B.Resize( A.Width(), A.Height() );
        Transpose( ASumFilt->LockedMatrix(), B.Matrix(), conjugate );
    }
}

template<typename T>
void TransposeContract
( const BlockMatrix<T>& A,
        BlockMatrix<T>& B, bool conjugate )
{
    DEBUG_CSE
    const Dist U = B.ColDist();
    const Dist V = B.RowDist();
    if( A.ColDist() == V && A.RowDist() == Partial(U) )
    {
        Transpose( A, B, conjugate );
    }
    else
    {
        unique_ptr<BlockMatrix<T>> 
            ASumFilt( B.ConstructTranspose(B.Grid(),B.Root()) );
        if( B.ColConstrained() )
            ASumFilt->AlignRowsWith( B, true );
        if( B.RowConstrained() )
            ASumFilt->AlignColsWith( B, true );
        Contract( A, *ASumFilt );
        if( !B.ColConstrained() )
            B.AlignColsWith( *ASumFilt, false );
        if( !B.RowConstrained() )
            B.AlignRowsWith( *ASumFilt, false );
        // We should have ensured that the alignments match
        B.Resize( A.Width(), A.Height() );
        Transpose( ASumFilt->LockedMatrix(), B.Matrix(), conjugate );
    }
}

template<typename T>
void AdjointContract( const ElementalMatrix<T>& A, ElementalMatrix<T>& B )
{
    DEBUG_CSE
    TransposeContract( A, B, true );
}

template<typename T>
void AdjointContract
( const BlockMatrix<T>& A,
        BlockMatrix<T>& B )
{
    DEBUG_CSE
    TransposeContract( A, B, true );
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void TransposeContract \
  ( const ElementalMatrix<T>& A, ElementalMatrix<T>& B, bool conjugate ); \
  EL_EXTERN template void TransposeContract \
  ( const BlockMatrix<T>& A, BlockMatrix<T>& B, bool conjugate ); \
  EL_EXTERN template void AdjointContract \
  ( const ElementalMatrix<T>& A, ElementalMatrix<T>& B ); \
  EL_EXTERN template void AdjointContract \
  ( const BlockMatrix<T>& A, BlockMatrix<T>& B );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_TRANSPOSECONTRACT_HPP
