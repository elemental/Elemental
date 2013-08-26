/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_LEGENDRE_HPP
#define ELEM_MATRICES_LEGENDRE_HPP

#include "elemental/matrices/Zeros.hpp"

namespace elem {

template<typename F> 
inline void
MakeLegendre( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("MakeLegendre");
#endif
    if( A.Height() != A.Width() )
        LogicError("Cannot make a non-square matrix Legendre");
    MakeZeros( A );

    const Int n = A.Width();
    for( Int j=0; j<n-1; ++j )
    {
        const F gamma = F(1) / Pow( F(2)*(j+1), F(2) );
        const F beta = F(1) / (2*Sqrt(F(1)-gamma));
        A.Set( j+1, j, beta );
        A.Set( j, j+1, beta );
    }
}

template<typename F,Distribution U,Distribution V>
inline void
MakeLegendre( DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    CallStackEntry cse("MakeLegendre");
#endif
    if( A.Height() != A.Width() )
        LogicError("Cannot make a non-square matrix Legendre");
    MakeZeros( A );

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
            if( j == i+1 || j == i-1 )
            {
                const Int k = Max( i, j );
                const F gamma = F(1) / Pow( F(2)*k, F(2) );
                const F beta = F(1) / (2*Sqrt(F(1)-gamma));
                A.SetLocal( iLoc, jLoc, beta );
            }
        }
    }
}

template<typename F> 
inline void
Legendre( Matrix<F>& A, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Legendre");
#endif
    A.ResizeTo( n, n );
    MakeLegendre( A );
}

template<typename F> 
inline Matrix<F>
Legendre( Int n )
{
    Matrix<F> A;
    Legendre( A, n );
    return A;
}

template<typename F,Distribution U,Distribution V> 
inline void
Legendre( DistMatrix<F,U,V>& A, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Legendre");
#endif
    A.ResizeTo( n, n );
    MakeLegendre( A );
}

template<typename F,Distribution U=MC,Distribution V=MR> 
inline DistMatrix<F,U,V>
Legendre( const Grid& g, Int n )
{
    DistMatrix<F,U,V> A(g);
    Legendre( A, n );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_LEGENDRE_HPP
