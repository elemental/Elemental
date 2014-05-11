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

#include ELEM_SETDIAGONAL_INC
#include ELEM_ZEROS_INC

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
    Zeros( A, n, n );
    SetDiagonal( A, 1,          1 );
    SetDiagonal( A, F(1)/F(4), -2 );
}

template<typename F,Dist U,Dist V>
inline void
Triangle( DistMatrix<F,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Triangle"))
    if( n < 3 )
        LogicError("Must be at least 3x3 to have a second-order symbol");
    Zeros( A, n, n );
    SetDiagonal( A, 1,          1 );
    SetDiagonal( A, F(1)/F(4), -2 );
}

template<typename F,Dist U,Dist V>
inline void
Triangle( BlockDistMatrix<F,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Triangle"))
    if( n < 3 )
        LogicError("Must be at least 3x3 to have a second-order symbol");
    Zeros( A, n, n );
    SetDiagonal( A, 1,          1 );
    SetDiagonal( A, F(1)/F(4), -2 );
}

} // namespace elem

#endif // ifndef ELEM_TRIANGLE_HPP
