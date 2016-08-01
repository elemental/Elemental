/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_TRANSPOSEAXPYCONTRACT_HPP
#define EL_BLAS_TRANSPOSEAXPYCONTRACT_HPP

namespace El {

template<typename T>
void TransposeAxpyContract
( T alpha, const ElementalMatrix<T>& A, 
                 ElementalMatrix<T>& B, bool conjugate )
{
    DEBUG_CSE
    const Dist U = B.ColDist();
    const Dist V = B.RowDist();
    if( A.ColDist() == V && A.RowDist() == U )
    {
        TransposeAxpy( alpha, A, B, conjugate );
    }
    else if( (A.ColDist() == V && A.RowDist() == Partial(U)) ||
             (A.ColDist() == V && A.RowDist() == Collect(U)) ||
             (A.RowDist() == U && A.ColDist() == Partial(V)) ||
             (A.RowDist() == U && A.ColDist() == Collect(V)) )
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

        // We should have ensured that the alignments are compatible
        TransposeAxpy( alpha, ASumFilt->LockedMatrix(), B.Matrix(), conjugate );
    }
    else
        LogicError("Incompatible distributions");
}

template<typename T>
void TransposeAxpyContract
( T alpha, const BlockMatrix<T>& A, 
                 BlockMatrix<T>& B, bool conjugate )
{
    DEBUG_CSE
    LogicError("Not yet implemented");
}

template<typename T>
void AdjointAxpyContract
( T alpha, const ElementalMatrix<T>& A,
                 ElementalMatrix<T>& B )
{
    DEBUG_CSE
    TransposeAxpyContract( alpha, A, B, true );
}

template<typename T>
void AdjointAxpyContract
( T alpha, const BlockMatrix<T>& A,
                 BlockMatrix<T>& B )
{
    DEBUG_CSE
    TransposeAxpyContract( alpha, A, B, true );
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void TransposeAxpyContract \
  ( T alpha, const ElementalMatrix<T>& A, \
                   ElementalMatrix<T>& B, bool conjugate ); \
  EL_EXTERN template void TransposeAxpyContract \
  ( T alpha, const BlockMatrix<T>& A, \
                   BlockMatrix<T>& B, bool conjugate ); \
  EL_EXTERN template void AdjointAxpyContract \
  ( T alpha, const ElementalMatrix<T>& A, ElementalMatrix<T>& B ); \
  EL_EXTERN template void AdjointAxpyContract \
  ( T alpha, const BlockMatrix<T>& A, BlockMatrix<T>& B );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_TRANSPOSEAXPYCONTRACT_HPP
