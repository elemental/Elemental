/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_FOURIER_HPP
#define ELEM_MATRICES_FOURIER_HPP

namespace elem {

template<typename R> 
inline void
MakeFourier( Matrix<Complex<R> >& A )
{
#ifndef RELEASE
    CallStackEntry cse("MakeFourier");
#endif
    typedef Complex<R> F;
    const Int m = A.Height();
    const Int n = A.Width();
    if( m != n )
        LogicError("Cannot make a non-square DFT matrix");

    const R pi = 4*Atan( R(1) );
    const R nSqrt = Sqrt( R(n) );
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            const R theta = -2*pi*i*j/n;
            const R realPart = cos(theta)/nSqrt;
            const R imagPart = sin(theta)/nSqrt;
            A.Set( i, j, F(realPart,imagPart) );
        }
    }
}

template<typename R,Distribution U,Distribution V>
inline void
MakeFourier( DistMatrix<Complex<R>,U,V>& A )
{
#ifndef RELEASE
    CallStackEntry cse("MakeFourier");
#endif
    typedef Complex<R> F;
    const Int m = A.Height();
    const Int n = A.Width();
    if( m != n )
        LogicError("Cannot make a non-square DFT matrix");

    const R pi = 4*Atan( R(1) );
    const R nSqrt = Sqrt( R(n) );
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
            const R theta = -2*pi*i*j/n;
            const R realPart = cos(theta)/nSqrt;
            const R imagPart = sin(theta)/nSqrt;
            A.SetLocal( iLoc, jLoc, F(realPart,imagPart) );
        }
    }
}

template<typename R>
inline void
Fourier( Matrix<Complex<R> >& A, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Fourier");
#endif
    A.ResizeTo( n, n );
    MakeFourier( A );
}

template<typename R>
inline Matrix<Complex<R> >
Fourier( Int n )
{
    Matrix<Complex<R> > A( n, n );
    MakeFourier( A );
    return A;
}

template<typename R,Distribution U,Distribution V>
inline void
Fourier( DistMatrix<Complex<R>,U,V>& A, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Fourier");
#endif
    A.ResizeTo( n, n );
    MakeFourier( A );
}

template<typename R,Distribution U=MC,Distribution V=MR>
inline DistMatrix<Complex<R>,U,V>
Fourier( const Grid& g, Int n )
{
    DistMatrix<Complex<R>,U,V> A( n, n, g );
    MakeFourier( A );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_FOURIER_HPP
