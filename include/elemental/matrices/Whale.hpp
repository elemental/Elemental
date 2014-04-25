/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_WHALE_HPP
#define ELEM_WHALE_HPP

#include "./Toeplitz.hpp"

namespace elem {

// The "whale matrix" is defined to have the symbol:
//   f(z) = -z^{-4} - (3+2i) z^{-3} + i z^{-2} + z^{-1} + 10 z + (3+i) z^2 +
//           4 z^3 + i z^4
// Please see 
//   A. Bottcher, "Infinite matrices and projection methods", 1996.

template<typename Real> 
inline void
Whale( Matrix<Complex<Real> >& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Whale"))
    if( n < 5 )
        LogicError("Must be at least 5x5 to have a fourth-order symbol");
    typedef Complex<Real> C;
    const Int numDiags = 2*n-1;
    std::vector<C> a( numDiags, 0 );
    const Int mainDiag = n-1;
    a[mainDiag-4] = -1;
    a[mainDiag-3] = C(-3,-2);
    a[mainDiag-2] = C( 0, 1);
    a[mainDiag-1] = 1;
    a[mainDiag-0] = 0;
    a[mainDiag+1] = 10;
    a[mainDiag+2] = C(3,1);
    a[mainDiag+3] = 4;
    a[mainDiag+4] = C(0,1);
    Toeplitz( A, n, n, a );
}

template<typename Real,Dist U,Dist V>
inline void
Whale( DistMatrix<Complex<Real>,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Whale"))
    if( n < 5 )
        LogicError("Must be at least 5x5 to have a fourth-order symbol");
    typedef Complex<Real> C;
    const Int numDiags = 2*n-1;
    std::vector<C> a( numDiags, 0 );
    const Int mainDiag = n-1;
    a[mainDiag-4] = -1;
    a[mainDiag-3] = C(-3,-2);
    a[mainDiag-2] = C( 0, 1);
    a[mainDiag-1] = 1;
    a[mainDiag-0] = 0;
    a[mainDiag+1] = 10;
    a[mainDiag+2] = C(3,1);
    a[mainDiag+3] = 4;
    a[mainDiag+4] = C(0,1);
    Toeplitz( A, n, n, a );
}

template<typename Real,Dist U,Dist V>
inline void
Whale( BlockDistMatrix<Complex<Real>,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Whale"))
    if( n < 5 )
        LogicError("Must be at least 5x5 to have a fourth-order symbol");
    typedef Complex<Real> C;
    const Int numDiags = 2*n-1;
    std::vector<C> a( numDiags, 0 );
    const Int mainDiag = n-1;
    a[mainDiag-4] = -1;
    a[mainDiag-3] = C(-3,-2);
    a[mainDiag-2] = C( 0, 1);
    a[mainDiag-1] = 1;
    a[mainDiag-0] = 0;
    a[mainDiag+1] = 10;
    a[mainDiag+2] = C(3,1);
    a[mainDiag+3] = 4;
    a[mainDiag+4] = C(0,1);
    Toeplitz( A, n, n, a );
}

#ifndef SWIG
template<typename Real> 
inline Matrix<Complex<Real>>
Whale( Int n )
{
    Matrix<Complex<Real>> A;
    Whale( A, n );
    return A;
}

template<typename Real,Dist U=MC,Dist V=MR>
inline DistMatrix<Complex<Real>,U,V>
Whale( const Grid& g, Int n )
{
    DistMatrix<Complex<Real>,U,V> A(g);
    Whale( A, n );
    return A;
}
#endif

} // namespace elem

#endif // ifndef ELEM_WHALE_HPP
