/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void Contract
( const AbstractDistMatrix<T>& A,
        AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Contract"))
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
        Zeros( B.Matrix(), B.LocalHeight(), B.LocalWidth() );
        AxpyContract( T(1), A, B );
    }
    else if( A.ColDist() == Partial(U) && A.RowDist() == V )
    {
        B.AlignAndResize
        ( A.ColAlign(), A.RowAlign(), A.Height(), A.Width(), false, false );
        Zeros( B.Matrix(), B.LocalHeight(), B.LocalWidth() );
        AxpyContract( T(1), A, B );
    }
    else if( A.ColDist() == U && A.RowDist() == Collect(V) )
    {
        B.AlignColsAndResize
        ( A.ColAlign(), A.Height(), A.Width(), false, false );
        Zeros( B.Matrix(), B.LocalHeight(), B.LocalWidth() );
        AxpyContract( T(1), A, B );
    }
    else if( A.ColDist() == Collect(U) && A.RowDist() == V )
    {
        B.AlignRowsAndResize
        ( A.RowAlign(), A.Height(), A.Width(), false, false );
        Zeros( B.Matrix(), B.LocalHeight(), B.LocalWidth() );
        AxpyContract( T(1), A, B );
    }
    else if( A.ColDist() == Collect(U) && A.RowDist() == Collect(V) )
    {
        Zeros( B, A.Height(), A.Width() );
        AxpyContract( T(1), A, B );
    }
    else
        LogicError("Incompatible distributions");
}

template<typename T>
void Contract
( const AbstractBlockDistMatrix<T>& A,
        AbstractBlockDistMatrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Contract"))
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
        Zeros( B.Matrix(), B.LocalHeight(), B.LocalWidth() );
        AxpyContract( T(1), A, B );
    }
    else if( A.ColDist() == Partial(U) && A.RowDist() == V )
    {
        B.AlignAndResize
        ( A.BlockHeight(), A.BlockWidth(),
          A.ColAlign(), A.RowAlign(), A.ColCut(), A.RowCut(),
          A.Height(), A.Width(), false, false );
        Zeros( B.Matrix(), B.LocalHeight(), B.LocalWidth() );
        AxpyContract( T(1), A, B );
    }
    else if( A.ColDist() == U && A.RowDist() == Collect(V) )
    {
        B.AlignColsAndResize
        ( A.BlockHeight(), A.ColAlign(), A.ColCut(), A.Height(), A.Width(),
          false, false );
        Zeros( B.Matrix(), B.LocalHeight(), B.LocalWidth() );
        AxpyContract( T(1), A, B );
    }
    else if( A.ColDist() == Collect(U) && A.RowDist() == V )
    {
        B.AlignRowsAndResize
        ( A.BlockWidth(), A.RowAlign(), A.RowCut(), A.Height(), A.Width(),
          false, false );
        Zeros( B.Matrix(), B.LocalHeight(), B.LocalWidth() );
        AxpyContract( T(1), A, B );
    }
    else if( A.ColDist() == Collect(U) && A.RowDist() == Collect(V) )
    {
        Zeros( B, A.Height(), A.Width() );
        AxpyContract( T(1), A, B );
    }
    else
        LogicError("Incompatible distributions");
}

#define PROTO(T) \
  template void Contract \
  ( const AbstractDistMatrix<T>& A, \
          AbstractDistMatrix<T>& B ); \
  template void Contract \
  ( const AbstractBlockDistMatrix<T>& A, \
          AbstractBlockDistMatrix<T>& B );

#include "El/macros/Instantiate.h"

} // namespace El
