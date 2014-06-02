/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_FOURIER_HPP
#define EL_FOURIER_HPP

namespace El {

template<typename Real> 
inline void
MakeFourier( Matrix<Complex<Real>>& A )
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

template<typename Real>
inline void
MakeFourier( AbstractDistMatrix<Complex<Real>>& A )
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
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            const Real theta = -2*pi*i*j/n;
            const Real realPart = Cos(theta)/nSqrt;
            const Real imagPart = Sin(theta)/nSqrt;
            A.SetLocal( iLoc, jLoc, Complex<Real>(realPart,imagPart) );
        }
    }
}

template<typename Real>
inline void
MakeFourier( AbstractBlockDistMatrix<Complex<Real>>& A )
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
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            const Real theta = -2*pi*i*j/n;
            const Real realPart = Cos(theta)/nSqrt;
            const Real imagPart = Sin(theta)/nSqrt;
            A.SetLocal( iLoc, jLoc, Complex<Real>(realPart,imagPart) );
        }
    }
}

template<typename Real>
inline void
Fourier( Matrix<Complex<Real>>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Fourier"))
    A.Resize( n, n );
    MakeFourier( A );
}

template<typename Real>
inline void
Fourier( AbstractDistMatrix<Complex<Real>>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Fourier"))
    A.Resize( n, n );
    MakeFourier( A );
}

template<typename Real>
inline void
Fourier( AbstractBlockDistMatrix<Complex<Real>>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Fourier"))
    A.Resize( n, n );
    MakeFourier( A );
}

} // namespace El

#endif // ifndef EL_FOURIER_HPP
