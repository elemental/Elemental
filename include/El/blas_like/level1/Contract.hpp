/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_CONTRACT_HPP
#define EL_BLAS_CONTRACT_HPP

namespace El {

template<typename T>
void Contract
( const ElementalMatrix<T>& A,
        ElementalMatrix<T>& B )
{
    EL_DEBUG_CSE
    AssertSameGrids( A, B );
    const Dist U = B.ColDist();
    const Dist V = B.RowDist();
    // TODO: Shorten this implementation?
    if( A.ColDist() == U && A.RowDist() == V )
    {
        Copy( A, B );
    }
    else if( A.ColDist() == U && A.RowDist() == Partial(V) )
    {
        B.AlignAndResize
        ( A.ColAlign(), A.RowAlign(), A.Height(), A.Width(), false, false );
        Zero( B.Matrix() );
        AxpyContract( T(1), A, B );
    }
    else if( A.ColDist() == Partial(U) && A.RowDist() == V )
    {
        B.AlignAndResize
        ( A.ColAlign(), A.RowAlign(), A.Height(), A.Width(), false, false );
        Zero( B.Matrix() );
        AxpyContract( T(1), A, B );
    }
    else if( A.ColDist() == U && A.RowDist() == Collect(V) )
    {
        B.AlignColsAndResize
        ( A.ColAlign(), A.Height(), A.Width(), false, false );
        Zero( B.Matrix() );
        AxpyContract( T(1), A, B );
    }
    else if( A.ColDist() == Collect(U) && A.RowDist() == V )
    {
        B.AlignRowsAndResize
        ( A.RowAlign(), A.Height(), A.Width(), false, false );
        Zero( B.Matrix() );
        AxpyContract( T(1), A, B );
    }
    else if( A.ColDist() == Collect(U) && A.RowDist() == Collect(V) )
    {
        B.Resize( A.Height(), A.Width() );
        Zero( B.Matrix() );
        AxpyContract( T(1), A, B );
    }
    else
        LogicError("Incompatible distributions");
}

template<typename T>
void Contract
( const BlockMatrix<T>& A,
        BlockMatrix<T>& B )
{
    EL_DEBUG_CSE
    AssertSameGrids( A, B );
    const Dist U = B.ColDist();
    const Dist V = B.RowDist();
    // TODO: Shorten this implementation?
    if( A.ColDist() == U && A.RowDist() == V )
    {
        Copy( A, B );
    }
    else if( A.ColDist() == U && A.RowDist() == Partial(V) )
    {
        B.AlignAndResize
        ( A.BlockHeight(), A.BlockWidth(),
          A.ColAlign(), A.RowAlign(), A.ColCut(), A.RowCut(),
          A.Height(), A.Width(), false, false );
        Zero( B.Matrix() );
        AxpyContract( T(1), A, B );
    }
    else if( A.ColDist() == Partial(U) && A.RowDist() == V )
    {
        B.AlignAndResize
        ( A.BlockHeight(), A.BlockWidth(),
          A.ColAlign(), A.RowAlign(), A.ColCut(), A.RowCut(),
          A.Height(), A.Width(), false, false );
        Zero( B.Matrix() );
        AxpyContract( T(1), A, B );
    }
    else if( A.ColDist() == U && A.RowDist() == Collect(V) )
    {
        B.AlignColsAndResize
        ( A.BlockHeight(), A.ColAlign(), A.ColCut(), A.Height(), A.Width(),
          false, false );
        Zero( B.Matrix() );
        AxpyContract( T(1), A, B );
    }
    else if( A.ColDist() == Collect(U) && A.RowDist() == V )
    {
        B.AlignRowsAndResize
        ( A.BlockWidth(), A.RowAlign(), A.RowCut(), A.Height(), A.Width(),
          false, false );
        Zero( B.Matrix() );
        AxpyContract( T(1), A, B );
    }
    else if( A.ColDist() == Collect(U) && A.RowDist() == Collect(V) )
    {
        B.Resize( A.Height(), A.Width() );
        Zero( B.Matrix() );
        AxpyContract( T(1), A, B );
    }
    else
        LogicError("Incompatible distributions");
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void Contract \
  ( const ElementalMatrix<T>& A, \
          ElementalMatrix<T>& B ); \
  EL_EXTERN template void Contract \
  ( const BlockMatrix<T>& A, \
          BlockMatrix<T>& B );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_CONTRACT_HPP
