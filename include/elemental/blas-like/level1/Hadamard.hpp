/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HADAMARD_HPP
#define ELEM_HADAMARD_HPP

// C(i,j) := A(i,j) B(i,j)

namespace elem {

template<typename T> 
inline void Hadamard( const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C )
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

template<typename T,Dist U,Dist V> 
inline void Hadamard
( const DistMatrix<T,U,V>& A, const DistMatrix<T,U,V>& B, DistMatrix<T,U,V>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Hadamard"))
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Hadamard product requires equal dimensions");
    if( A.Grid() != B.Grid() )
        LogicError("A and B must have the same grids");
    if( A.ColAlign() != B.ColAlign() || A.RowAlign() != B.RowAlign() )
        LogicError("A and B must be aligned");
    C.AlignWith( A );
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

} // namespace elem

#endif // ifndef ELEM_HADAMARD_HPP
