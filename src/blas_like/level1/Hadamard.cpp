/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// C(i,j) := A(i,j) B(i,j)

namespace El {

template<typename T> 
void Hadamard( const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Hadamard"))
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Hadamard product requires equal dimensions");
    C.Resize( A.Height(), A.Width() );

    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            C.Set( i, j, A.Get(i,j)*B.Get(i,j) );
}

template<typename T> 
void Hadamard
( const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B, 
        AbstractDistMatrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Hadamard"))
    const DistData ADistData = A.DistData();
    const DistData BDistData = B.DistData();
    DistData CDistData = C.DistData();
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Hadamard product requires equal dimensions");
    AssertSameGrids( A, B );
    if( ADistData.colDist != BDistData.colDist ||
        ADistData.rowDist != BDistData.rowDist ||
        BDistData.colDist != CDistData.colDist ||
        BDistData.rowDist != CDistData.rowDist )
        LogicError("A, B, and C must share the same distribution");
    if( A.ColAlign() != B.ColAlign() || A.RowAlign() != B.RowAlign() )
        LogicError("A and B must be aligned");
    C.AlignWith( A.DistData() );
    C.Resize( A.Height(), A.Width() );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const T alpha = A.GetLocal(iLoc,jLoc); 
            const T beta = B.GetLocal(iLoc,jLoc);
            C.SetLocal( iLoc, jLoc, alpha*beta );
        }
    }
}

#define PROTO(T) \
  template void Hadamard \
  ( const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C ); \
  template void Hadamard \
  ( const AbstractDistMatrix<T>& A, \
    const AbstractDistMatrix<T>& B, \
          AbstractDistMatrix<T>& C );

#include "El/macros/Instantiate.h"

} // namespace El
