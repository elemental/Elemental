/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_WHALE_HPP
#define EL_WHALE_HPP



namespace El {

// The "whale matrix" is defined to have the symbol:
//   f(z) = -z^{-4} - (3+2i) z^{-3} + i z^{-2} + z^{-1} + 10 z + (3+i) z^2 +
//           4 z^3 + i z^4
// Please see 
//   A. Bottcher, "Infinite matrices and projection methods", 1996.

template<typename Real> 
inline void
Whale( Matrix<Complex<Real>>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Whale"))
    if( n < 5 )
        LogicError("Must be at least 5x5 to have a fourth-order symbol");
    typedef Complex<Real> C;
    Zeros( A, n, n );
    SetDiagonal( A, -1,        4 );
    SetDiagonal( A, C(-3,-2),  3 );
    SetDiagonal( A, C( 0, 1),  2 );
    SetDiagonal( A,  1,        1 );
    SetDiagonal( A, 10,       -1 );
    SetDiagonal( A, C( 3, 1), -2 );
    SetDiagonal( A,  4,       -3 );
    SetDiagonal( A, C( 0, 1), -4 );
}

template<typename Real>
inline void
Whale( AbstractDistMatrix<Complex<Real>>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Whale"))
    if( n < 5 )
        LogicError("Must be at least 5x5 to have a fourth-order symbol");
    typedef Complex<Real> C;
    Zeros( A, n, n );
    SetDiagonal( A, -1,        4 );
    SetDiagonal( A, C(-3,-2),  3 );
    SetDiagonal( A, C( 0, 1),  2 );
    SetDiagonal( A,  1,        1 );
    SetDiagonal( A, 10,       -1 );
    SetDiagonal( A, C( 3, 1), -2 );
    SetDiagonal( A,  4,       -3 );
    SetDiagonal( A, C( 0, 1), -4 );
}

template<typename Real>
inline void
Whale( AbstractBlockDistMatrix<Complex<Real>>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Whale"))
    if( n < 5 )
        LogicError("Must be at least 5x5 to have a fourth-order symbol");
    typedef Complex<Real> C;
    Zeros( A, n, n );
    SetDiagonal( A, -1,        4 );
    SetDiagonal( A, C(-3,-2),  3 );
    SetDiagonal( A, C( 0, 1),  2 );
    SetDiagonal( A,  1,        1 );
    SetDiagonal( A, 10,       -1 );
    SetDiagonal( A, C( 3, 1), -2 );
    SetDiagonal( A,  4,       -3 );
    SetDiagonal( A, C( 0, 1), -4 );
}

} // namespace El

#endif // ifndef EL_WHALE_HPP
