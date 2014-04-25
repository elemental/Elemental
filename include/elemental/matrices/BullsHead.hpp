/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BULLSHEAD_HPP
#define ELEM_BULLSHEAD_HPP

#include "./Toeplitz.hpp"

namespace elem {

// The "bull's head matrix" is defined to have the symbol:
//   f(z) = 2i z^{-1} + z^2 + 7/10 z^3.
// Please see 
//   L. Reichel and L. N. Trefethen, "Eigenvalues and pseudo-eigenvalues of 
//   "Toeplitz matrices", Linear Algebra Appl., 1992.

template<typename Real> 
inline void
BullsHead( Matrix<Complex<Real> >& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("BullsHead"))
    if( n < 4 )
        LogicError("Must be at least 4x4 to have a third-order symbol");
    typedef Complex<Real> C;
    const Int numDiags = 2*n-1;
    std::vector<C> a( numDiags, 0 );
    const Int mainDiag = n-1;
    a[mainDiag-1] = C(0,2);
    a[mainDiag-0] = 0;
    a[mainDiag+1] = 0;
    a[mainDiag+2] = 1;
    a[mainDiag+3] = Real(7)/Real(10);
    Toeplitz( A, n, n, a );
}

template<typename Real,Dist U,Dist V>
inline void
BullsHead( DistMatrix<Complex<Real>,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("BullsHead"))
    if( n < 4 )
        LogicError("Must be at least 4x4 to have a third-order symbol");
    typedef Complex<Real> C;
    const Int numDiags = 2*n-1;
    std::vector<C> a( numDiags, 0 );
    const Int mainDiag = n-1;
    a[mainDiag-1] = C(0,2);
    a[mainDiag-0] = 0;
    a[mainDiag+1] = 0;
    a[mainDiag+2] = 1;
    a[mainDiag+3] = Real(7)/Real(10);
    Toeplitz( A, n, n, a );
}

template<typename Real,Dist U,Dist V>
inline void
BullsHead( BlockDistMatrix<Complex<Real>,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("BullsHead"))
    if( n < 4 )
        LogicError("Must be at least 4x4 to have a third-order symbol");
    typedef Complex<Real> C;
    const Int numDiags = 2*n-1;
    std::vector<C> a( numDiags, 0 );
    const Int mainDiag = n-1;
    a[mainDiag-1] = C(0,2);
    a[mainDiag-0] = 0;
    a[mainDiag+1] = 0;
    a[mainDiag+2] = 1;
    a[mainDiag+3] = Real(7)/Real(10);
    Toeplitz( A, n, n, a );
}

#ifndef SWIG
template<typename Real> 
inline Matrix<Complex<Real>>
BullsHead( Int n )
{
    Matrix<Complex<Real>> A;
    BullsHead( A, n );
    return A;
}

template<typename Real,Dist U=MC,Dist V=MR>
inline DistMatrix<Complex<Real>,U,V>
BullsHead( const Grid& g, Int n )
{
    DistMatrix<Complex<Real>,U,V> A(g);
    BullsHead( A, n );
    return A;
}

template<typename Real,Dist U=MC,Dist V=MR>
inline BlockDistMatrix<Complex<Real>,U,V>
BullsHead( const Grid& g, Int n )
{
    BlockDistMatrix<Complex<Real>,U,V> A(g);
    BullsHead( A, n );
    return A;
}
#endif

} // namespace elem

#endif // ifndef ELEM_BULLSHEAD_HPP
