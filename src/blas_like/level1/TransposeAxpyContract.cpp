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
void TransposeAxpyContract
( T alpha, const AbstractDistMatrix<T>& A, 
                 AbstractDistMatrix<T>& B, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("TransposeAxpyContract"))
    const Dist U = B.ColDist();
    const Dist V = B.RowDist();
    if( A.ColDist() == V && A.RowDist() == U )
    {
        TransposeAxpy( alpha, A, B, conjugate );
    }
    else if( A.ColDist() == V && A.RowDist() == Partial(U) )
    {
        std::unique_ptr<AbstractDistMatrix<T>>
          ASumFilt( B.ConstructTranspose(B.Grid(),B.Root()) );
        if( B.ColConstrained() )
            ASumFilt->AlignRowsWith( B, false );
        if( B.RowConstrained() )
            ASumFilt->AlignColsWith( B, false );
        Contract( A, *ASumFilt );
        if( !B.ColConstrained() )
            B.AlignColsWith( *ASumFilt, false );
        if( !B.RowConstrained() )
            B.AlignRowsWith( *ASumFilt, false );
        TransposeAxpy( alpha, ASumFilt->LockedMatrix(), B.Matrix(), conjugate );
    }
    else if( A.ColDist() == V && A.RowDist() == Collect(U) )
    {
        std::unique_ptr<AbstractDistMatrix<T>>
          ASumFilt( B.ConstructTranspose(B.Grid(),B.Root()) );
        if( B.ColConstrained() )
            ASumFilt->AlignRowsWith( B, false );
        if( B.RowConstrained() )
            ASumFilt->AlignColsWith( B, false );
        Contract( A, *ASumFilt );
        if( !B.ColConstrained() )
            B.AlignColsWith( *ASumFilt, false );
        if( !B.RowConstrained() )
            B.AlignRowsWith( *ASumFilt, false );
        TransposeAxpy( alpha, ASumFilt->LockedMatrix(), B.Matrix(), conjugate );
    }
    else
        LogicError("Incompatible distributions");
}

template<typename T>
void TransposeAxpyContract
( T alpha, const AbstractBlockDistMatrix<T>& A, 
                 AbstractBlockDistMatrix<T>& B, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("TransposeAxpyContract"))
    LogicError("Not yet implemented");
}

#define PROTO(T) \
  template void TransposeAxpyContract \
  ( T alpha, const AbstractDistMatrix<T>& A, \
                   AbstractDistMatrix<T>& B, bool conjugate ); \
  template void TransposeAxpyContract \
  ( T alpha, const AbstractBlockDistMatrix<T>& A, \
                   AbstractBlockDistMatrix<T>& B, bool conjugate );

#include "El/macros/Instantiate.h"

} // namespace El
