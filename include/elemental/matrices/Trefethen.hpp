/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TREFETHEN_HPP
#define ELEM_TREFETHEN_HPP

#include "./Toeplitz.hpp"

namespace elem {

// The symbol of the matrix described at the beginning of Chapter II of 
// Lloyd N. Trefethen and Mark Embree's "Spectra and Pseudospectra" is
//     f(z) = 2 z^{-3} - z^{-2} + 2i z^{-1} - 4 z^2 - 2i z^3.
// For lack of a better name, we will refer to such a Toeplitz matrix as a 
// Trefethen matrix.

template<typename Real> 
inline void
Trefethen( Matrix<Complex<Real> >& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Trefethen"))
    if( n < 4 )
        LogicError("Must be at least 4x4 to have a third-order symbol");
    typedef Complex<Real> C;
    const Int numDiags = 2*n-1;
    std::vector<C> a( numDiags, 0 );
    const Int mainDiag = n-1;
    a[mainDiag-3] = 2;
    a[mainDiag-2] = -1;
    a[mainDiag-1] = C(0,2);
    a[mainDiag-0] = 0;
    a[mainDiag+1] = 0;
    a[mainDiag+2] = -4;
    a[mainDiag+3] = C(0,-2);
    Toeplitz( A, n, n, a );
}

template<typename Real,Dist U,Dist V>
inline void
Trefethen( DistMatrix<Complex<Real>,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Trefethen"))
    if( n < 4 )
        LogicError("Must be at least 4x4 to have a third-order symbol");
    typedef Complex<Real> C;
    const Int numDiags = 2*n-1;
    std::vector<C> a( numDiags, 0 );
    const Int mainDiag = n-1;
    a[mainDiag-3] = 2;
    a[mainDiag-2] = -1;
    a[mainDiag-1] = C(0,2);
    a[mainDiag-0] = 0;
    a[mainDiag+1] = 0;
    a[mainDiag+2] = -4;
    a[mainDiag+3] = C(0,-2);
    Toeplitz( A, n, n, a );
}

template<typename Real,Dist U,Dist V>
inline void
Trefethen( BlockDistMatrix<Complex<Real>,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Trefethen"))
    if( n < 4 )
        LogicError("Must be at least 4x4 to have a third-order symbol");
    typedef Complex<Real> C;
    const Int numDiags = 2*n-1;
    std::vector<C> a( numDiags, 0 );
    const Int mainDiag = n-1;
    a[mainDiag-3] = 2;
    a[mainDiag-2] = -1;
    a[mainDiag-1] = C(0,2);
    a[mainDiag-0] = 0;
    a[mainDiag+1] = 0;
    a[mainDiag+2] = -4;
    a[mainDiag+3] = C(0,-2);
    Toeplitz( A, n, n, a );
}

#ifndef SWIG
template<typename Real> 
inline Matrix<Complex<Real>>
Trefethen( Int n )
{
    Matrix<Complex<Real>> A;
    Trefethen( A, n );
    return A;
}

template<typename Real,Dist U=MC,Dist V=MR>
inline DistMatrix<Complex<Real>,U,V>
Trefethen( const Grid& g, Int n )
{
    DistMatrix<Complex<Real>,U,V> A(g);
    Trefethen( A, n );
    return A;
}
#endif

} // namespace elem

#endif // ifndef ELEM_TREFETHEN_HPP
