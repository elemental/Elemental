/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_HILBERT_HPP
#define ELEM_MATRICES_HILBERT_HPP

namespace elem {

template<typename F> 
inline void
MakeHilbert( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("MakeHilbert");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    if( m != n )
        LogicError("Cannot make a non-square matrix Hilbert");

    const F one = F(1);
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, one/F(i+j+1) );
}

template<typename F,Distribution U,Distribution V>
inline void
MakeHilbert( DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    CallStackEntry cse("MakeHilbert");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    if( m != n )
        LogicError("Cannot make a non-square matrix Hilbert");

    const F one = F(1);
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int colShift = A.ColShift();
    const Int rowShift = A.RowShift();
    const Int colStride = A.ColStride();
    const Int rowStride = A.RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*colStride;
            A.SetLocal( iLoc, jLoc, one/F(i+j+1) );
        }
    }
}

template<typename F>
inline void
Hilbert( Matrix<F>& A, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Hilbert");
#endif
    A.ResizeTo( n, n );
    MakeHilbert( A );
}

template<typename F>
inline Matrix<F>
Hilbert( Int n )
{
    Matrix<F> A( n, n );
    MakeHilbert( A );
    return A;
}

template<typename F,Distribution U,Distribution V>
inline void
Hilbert( DistMatrix<F,U,V>& A, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Hilbert");
#endif
    A.ResizeTo( n, n );
    MakeHilbert( A );
}

template<typename F,Distribution U=MC,Distribution V=MR>
inline DistMatrix<F,U,V>
Hilbert( const Grid& g, Int n )
{
    DistMatrix<F,U,V> A( n, n, g );
    MakeHilbert( A );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_HILBERT_HPP
