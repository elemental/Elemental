/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace transpose {

// (Partial(V),U) |-> (U,V)
template<typename T>
void PartialRowFilter
( const AbstractDistMatrix<T>& A, 
        AbstractDistMatrix<T>& B, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("transpose::PartialRowFilter");
        if( A.ColDist() != Partial(B.RowDist()) ||
            A.RowDist() != B.ColDist() )
            LogicError("Incompatible distributions");
    )
    std::unique_ptr<AbstractDistMatrix<T>>
        AFilt( B.ConstructTranspose(B.Grid(),B.Root()) );
    if( B.ColConstrained() )
        AFilt->AlignRowsWith( B, false );
    if( B.RowConstrained() )
        AFilt->AlignColsWith( B, false );
    Copy( A, *AFilt );
    if( !B.ColConstrained() )
        B.AlignColsWith( *AFilt, false );
    if( !B.RowConstrained() )
        B.AlignRowsWith( *AFilt, false );
    B.Resize( A.Width(), A.Height() );
    Transpose( AFilt->LockedMatrix(), B.Matrix(), conjugate );
}

template<typename T>
void PartialRowFilter
( const AbstractBlockDistMatrix<T>& A, 
        AbstractBlockDistMatrix<T>& B, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("transpose::PartialRowFilter");
        if( A.ColDist() != Partial(B.RowDist()) ||
            A.RowDist() != B.ColDist() )
            LogicError("Incompatible distributions");
    )
    std::unique_ptr<AbstractBlockDistMatrix<T>>
        AFilt( B.ConstructTranspose(B.Grid(),B.Root()) );
    if( B.ColConstrained() )
        AFilt->AlignRowsWith( B, false );
    if( B.RowConstrained() )
        AFilt->AlignColsWith( B, false );
    Copy( A, *AFilt );
    if( !B.ColConstrained() )
        B.AlignColsWith( *AFilt, false );
    if( !B.RowConstrained() )
        B.AlignRowsWith( *AFilt, false );
    B.Resize( A.Width(), A.Height() );
    Transpose( AFilt->LockedMatrix(), B.Matrix(), conjugate );
}

#define PROTO(T) \
  template void PartialRowFilter \
  ( const AbstractDistMatrix<T>& A, \
          AbstractDistMatrix<T>& B, bool conjugate ); \
  template void PartialRowFilter \
  ( const AbstractBlockDistMatrix<T>& A, \
          AbstractBlockDistMatrix<T>& B, bool conjugate );

#include "El/macros/Instantiate.h"

} // namespace transpose
} // namespace El
