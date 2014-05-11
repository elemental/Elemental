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

#include ELEM_SETDIAGONAL_INC
#include ELEM_ZEROS_INC

namespace elem {

// The "bull's head matrix" is defined to have the symbol:
//   f(z) = 2i z^{-1} + z^2 + 7/10 z^3.
// Please see 
//   L. Reichel and L. N. Trefethen, "Eigenvalues and pseudo-eigenvalues of 
//   "Toeplitz matrices", Linear Algebra Appl., 1992.

template<typename Real> 
inline void
BullsHead( Matrix<Complex<Real>>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("BullsHead"))
    if( n < 4 )
        LogicError("Must be at least 4x4 to have a third-order symbol");
    typedef Complex<Real> C;
    Zeros( A, n, n );
    SetDiagonal( A, C(0,2),            1 );
    SetDiagonal( A, 1,                -2 );
    SetDiagonal( A, Real(7)/Real(10), -3 );
}

template<typename Real,Dist U,Dist V>
inline void
BullsHead( DistMatrix<Complex<Real>,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("BullsHead"))
    if( n < 4 )
        LogicError("Must be at least 4x4 to have a third-order symbol");
    typedef Complex<Real> C;
    Zeros( A, n, n );
    SetDiagonal( A, C(0,2),            1 );
    SetDiagonal( A, 1,                -2 );
    SetDiagonal( A, Real(7)/Real(10), -3 );
}

template<typename Real,Dist U,Dist V>
inline void
BullsHead( BlockDistMatrix<Complex<Real>,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("BullsHead"))
    if( n < 4 )
        LogicError("Must be at least 4x4 to have a third-order symbol");
    typedef Complex<Real> C;
    Zeros( A, n, n );
    SetDiagonal( A, C(0,2),            1 );
    SetDiagonal( A, 1,                -2 );
    SetDiagonal( A, Real(7)/Real(10), -3 );
}

} // namespace elem

#endif // ifndef ELEM_BULLSHEAD_HPP
