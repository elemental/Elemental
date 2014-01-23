/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_FOURIER_HPP
#define ELEM_MATRICES_FOURIER_HPP

namespace elem {

template<typename Real> 
inline void
MakeFourier( Matrix<Complex<Real> >& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeFourier"))
    const Int m = A.Height();
    const Int n = A.Width();
    if( m != n )
        LogicError("Cannot make a non-square DFT matrix");

    const Real pi = 4*Atan( Real(1) );
    const Real nSqrt = Sqrt( Real(n) );
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            const Real theta = -2*pi*i*j/n;
            const Real realPart = Cos(theta)/nSqrt;
            const Real imagPart = Sin(theta)/nSqrt;
            A.Set( i, j, Complex<Real>(realPart,imagPart) );
        }
    }
}

template<typename Real,Distribution U,Distribution V>
inline void
MakeFourier( DistMatrix<Complex<Real>,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeFourier"))
    const Int m = A.Height();
    const Int n = A.Width();
    if( m != n )
        LogicError("Cannot make a non-square DFT matrix");

    const Real pi = 4*Atan( Real(1) );
    const Real nSqrt = Sqrt( Real(n) );
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
            const Real theta = -2*pi*i*j/n;
            const Real realPart = Cos(theta)/nSqrt;
            const Real imagPart = Sin(theta)/nSqrt;
            A.SetLocal( iLoc, jLoc, Complex<Real>(realPart,imagPart) );
        }
    }
}

template<typename Real>
inline void
Fourier( Matrix<Complex<Real> >& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Fourier"))
    A.Resize( n, n );
    MakeFourier( A );
}

template<typename Real>
inline Matrix<Complex<Real> >
Fourier( Int n )
{
    Matrix<Complex<Real>> A( n, n );
    MakeFourier( A );
    return A;
}

template<typename Real,Distribution U,Distribution V>
inline void
Fourier( DistMatrix<Complex<Real>,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Fourier"))
    A.Resize( n, n );
    MakeFourier( A );
}

template<typename Real,Distribution U=MC,Distribution V=MR>
inline DistMatrix<Complex<Real>,U,V>
Fourier( const Grid& g, Int n )
{
    DistMatrix<Complex<Real>,U,V> A( n, n, g );
    MakeFourier( A );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_FOURIER_HPP
