/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TRIANGLE_HPP
#define ELEM_TRIANGLE_HPP

#include "./Toeplitz.hpp"

namespace elem {

// The "triangle matrix" is defined to have the symbol:
//   f(z) = z^{-1} + 1/4 z^2.
// Please see 
//   L. Reichel and L. N. Trefethen, "Eigenvalues and pseudo-eigenvalues of 
//   "Toeplitz matrices", Linear Algebra Appl., 1992.

template<typename F> 
inline void
Triangle( Matrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Triangle"))
    if( n < 3 )
        LogicError("Must be at least 3x3 to have a second-order symbol");
    const Int numDiags = 2*n-1;
    std::vector<F> a( numDiags, 0 );
    a[n-2] = 1;
    a[n-1] = 0;
    a[n  ] = 0;
    a[n+1] = F(1)/F(4);
    Toeplitz( A, n, n, a );
}

template<typename F,Dist U,Dist V>
inline void
Triangle( DistMatrix<F,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Triangle"))
    if( n < 3 )
        LogicError("Must be at least 3x3 to have a second-order symbol");
    const Int numDiags = 2*n-1;
    std::vector<F> a( numDiags, 0 );
    a[n-2] = 1;
    a[n-1] = 0;
    a[n  ] = 0;
    a[n+1] = F(1)/F(4);
    Toeplitz( A, n, n, a );
}

template<typename F,Dist U,Dist V>
inline void
Triangle( BlockDistMatrix<F,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Triangle"))
    if( n < 3 )
        LogicError("Must be at least 3x3 to have a second-order symbol");
    const Int numDiags = 2*n-1;
    std::vector<F> a( numDiags, 0 );
    a[n-2] = 1;
    a[n-1] = 0;
    a[n  ] = 0;
    a[n+1] = F(1)/F(4);
    Toeplitz( A, n, n, a );
}

#ifndef SWIG
template<typename F> 
inline Matrix<F>
Triangle( Int n )
{
    Matrix<F> A;
    Triangle( A, n );
    return A;
}

template<typename F,Dist U=MC,Dist V=MR>
inline DistMatrix<F,U,V>
Triangle( const Grid& g, Int n )
{
    DistMatrix<F,U,V> A(g);
    Triangle( A, n );
    return A;
}
#endif

} // namespace elem

#endif // ifndef ELEM_TRIANGLE_HPP
