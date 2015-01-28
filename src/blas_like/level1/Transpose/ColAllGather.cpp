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

// (U,V) |-> (V,Collect(U))
template<typename T>
void ColAllGather
( const AbstractDistMatrix<T>& A, 
        AbstractDistMatrix<T>& B, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("transpose::ColAllGather");
        if( B.ColDist() != A.RowDist() ||
            B.RowDist() != Collect(A.ColDist()) )
            LogicError("Incompatible distributions");
    )
    std::unique_ptr<AbstractDistMatrix<T>> ATrans
    ( A.ConstructTranspose(A.Grid(),A.Root()) );
    ATrans->AlignWith( A );
    ATrans->Resize( A.Width(), A.Height() );
    Transpose( A.LockedMatrix(), ATrans->Matrix(), conjugate );
    Copy( *ATrans, B );
}

template<typename T>
void ColAllGather
( const AbstractBlockDistMatrix<T>& A, 
        AbstractBlockDistMatrix<T>& B, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("transpose::ColAllGather");
        if( B.ColDist() != A.RowDist() ||
            B.RowDist() != Collect(A.ColDist()) )
            LogicError("Incompatible distributions");
    )
    std::unique_ptr<AbstractBlockDistMatrix<T>> ATrans
    ( A.ConstructTranspose(A.Grid(),A.Root()) );
    ATrans->AlignWith( A );
    ATrans->Resize( A.Width(), A.Height() );
    Transpose( A.LockedMatrix(), ATrans->Matrix(), conjugate );
    Copy( *ATrans, B );
}

#define PROTO(T) \
  template void ColAllGather \
  ( const AbstractDistMatrix<T>& A, \
          AbstractDistMatrix<T>& B, bool conjugate ); \
  template void ColAllGather \
  ( const AbstractBlockDistMatrix<T>& A, \
          AbstractBlockDistMatrix<T>& B, bool conjugate );

#include "El/macros/Instantiate.h"

} // namespace transpose
} // namespace El
