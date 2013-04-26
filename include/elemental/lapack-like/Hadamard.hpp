/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HADAMARD_HPP
#define LAPACK_HADAMARD_HPP

//
// C(i,j) := A(i,j) B(i,j)
//

namespace elem {

template<typename T> 
inline void Hadamard( const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C )
{
#ifndef RELEASE
    CallStackEntry entry("Hadamard");
#endif
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        throw std::logic_error("Hadamard product requires equal dimensions");
    C.ResizeTo( A.Height(), A.Width() );

    const int height = A.Height();
    const int width = A.Width();
    for( int j=0; j<width; ++j )
        for( int i=0; i<height; ++i )
            C.Set( i, j, A.Get(i,j)*B.Get(i,j) );
}

template<typename T,Distribution U,Distribution V> 
inline void Hadamard
( const DistMatrix<T,U,V>& A, const DistMatrix<T,U,V>& B, DistMatrix<T,U,V>& C )
{
#ifndef RELEASE
    CallStackEntry entry("Hadamard");
#endif
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        throw std::logic_error("Hadamard product requires equal dimensions");
    if( A.Grid() != B.Grid() )
        throw std::logic_error("A and B must have the same grids");
    if( A.ColAlignment() != B.ColAlignment() || 
        A.RowAlignment() != B.RowAlignment() )
        throw std::logic_error("A and B must be aligned");
    const Grid& g = A.Grid();
    C.SetGrid( g );
    C.AlignWith( A );
    C.ResizeTo( A.Height(), A.Width() );

    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const T alpha = A.GetLocal(iLocal,jLocal); 
            const T beta = B.GetLocal(iLocal,jLocal);
            C.SetLocal( iLocal, jLocal, alpha*beta );
        }
    }
}

} // namespace elem

#endif // ifndef LAPACK_HADAMARD_HPP
