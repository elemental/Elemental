/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_FIEDLER_HPP
#define ELEM_MATRICES_FIEDLER_HPP

namespace elem {

template<typename F> 
inline void
Fiedler( Matrix<F>& A, const std::vector<F>& c )
{
#ifndef RELEASE
    CallStackEntry cse("Fiedler");
#endif
    const Int n = c.size();
    A.ResizeTo( n, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<n; ++i )
            A.Set( i, j, Abs(c[i]-c[j]) );
}

template<typename F> 
inline Matrix<F>
Fiedler( const std::vector<F>& c )
{
    Matrix<F> A;
    Fiedler( A, c ); 
    return A;
}

template<typename F,Distribution U,Distribution V>
inline void
Fiedler( DistMatrix<F,U,V>& A, const std::vector<F>& c )
{
#ifndef RELEASE
    CallStackEntry cse("Fiedler");
#endif
    const Int n = c.size();
    A.ResizeTo( n, n );
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
            A.SetLocal( iLoc, jLoc, Abs(c[i]-c[j]) );
        }
    }
}

template<typename F,Distribution U=MC,Distribution V=MR>
inline DistMatrix<F,U,V>
Fiedler( const Grid& g, const std::vector<F>& c )
{
    DistMatrix<F,U,V> A(g);
    Fiedler( A, c );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_FIEDLER_HPP
