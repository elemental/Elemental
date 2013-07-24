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
Fourier( Matrix<Complex<R> >& A, int n )
{
#ifndef RELEASE
    CallStackEntry entry("Fourier");
#endif
    A.ResizeTo( n, n );
    MakeFourier( A );
}

template<typename R,Distribution U,Distribution V>
inline void
Fourier( DistMatrix<Complex<R>,U,V>& A, int n )
{
#ifndef RELEASE
    CallStackEntry entry("Fourier");
#endif
    A.ResizeTo( n, n );
    MakeFourier( A );
}

template<typename R> 
inline void
MakeFourier( Matrix<Complex<R> >& A )
{
#ifndef RELEASE
    CallStackEntry entry("MakeFourier");
#endif
    typedef Complex<R> F;
    const int m = A.Height();
    const int n = A.Width();
    if( m != n )
        throw std::logic_error("Cannot make a non-square DFT matrix");

    const R pi = 4*Atan( R(1) );
    const R nSqrt = Sqrt( R(n) );
    for( int j=0; j<n; ++j )
    {
        for( int i=0; i<m; ++i )
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
    CallStackEntry entry("MakeFourier");
#endif
    typedef Complex<R> F;
    const int m = A.Height();
    const int n = A.Width();
    if( m != n )
        throw std::logic_error("Cannot make a non-square DFT matrix");

    const R pi = 4*Atan( R(1) );
    const R nSqrt = Sqrt( R(n) );
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();
    const int colStride = A.ColStride();
    const int rowStride = A.RowStride();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const int i = colShift + iLoc*colStride;
            const R theta = -2*pi*i*j/n;
            const R realPart = cos(theta)/nSqrt;
            const R imagPart = sin(theta)/nSqrt;
            A.SetLocal( iLoc, jLoc, F(realPart,imagPart) );
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_FOURIER_HPP
