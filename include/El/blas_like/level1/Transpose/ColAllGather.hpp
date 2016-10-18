/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_TRANSPOSE_COLALLGATHER_HPP
#define EL_BLAS_TRANSPOSE_COLALLGATHER_HPP

namespace El {
namespace transpose {

// (U,V) |-> (V,Collect(U))
template<typename T>
void ColAllGather
( const ElementalMatrix<T>& A, 
        ElementalMatrix<T>& B, bool conjugate )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( B.ColDist() != A.RowDist() ||
          B.RowDist() != Collect(A.ColDist()) )
          LogicError("Incompatible distributions");
    )
    unique_ptr<ElementalMatrix<T>> ATrans
    ( A.ConstructTranspose(A.Grid(),A.Root()) );
    ATrans->AlignWith( A );
    ATrans->Resize( A.Width(), A.Height() );
    Transpose( A.LockedMatrix(), ATrans->Matrix(), conjugate );
    Copy( *ATrans, B );
}

template<typename T>
void ColAllGather
( const BlockMatrix<T>& A, 
        BlockMatrix<T>& B, bool conjugate )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( B.ColDist() != A.RowDist() ||
          B.RowDist() != Collect(A.ColDist()) )
          LogicError("Incompatible distributions");
    )
    unique_ptr<BlockMatrix<T>> ATrans
    ( A.ConstructTranspose(A.Grid(),A.Root()) );
    ATrans->AlignWith( A );
    ATrans->Resize( A.Width(), A.Height() );
    Transpose( A.LockedMatrix(), ATrans->Matrix(), conjugate );
    Copy( *ATrans, B );
}

} // namespace transpose
} // namespace El

#endif // ifndef EL_BLAS_TRANSPOSE_COLALLGATHER_HPP
