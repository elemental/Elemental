/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

// The "whale matrix" is defined to have the symbol:
//   f(z) = -z^{-4} - (3+2i) z^{-3} + i z^{-2} + z^{-1} + 10 z + (3+i) z^2 +
//           4 z^3 + i z^4
// Please see 
//   A. Bottcher, "Infinite matrices and projection methods", 1996.

template<typename Real> 
void Whale( Matrix<Complex<Real>>& A, Int n )
{
    EL_DEBUG_CSE
    if( n < 5 )
        LogicError("Must be at least 5x5 to have a fourth-order symbol");
    typedef Complex<Real> C;
    Zeros( A, n, n );
    FillDiagonal( A, C(-1),     4 );
    FillDiagonal( A, C(-3,-2),  3 );
    FillDiagonal( A, C( 0, 1),  2 );
    FillDiagonal( A, C(1),      1 );
    FillDiagonal( A, C(10),    -1 );
    FillDiagonal( A, C( 3, 1), -2 );
    FillDiagonal( A, C(4),     -3 );
    FillDiagonal( A, C( 0, 1), -4 );
}

template<typename Real>
void Whale( AbstractDistMatrix<Complex<Real>>& A, Int n )
{
    EL_DEBUG_CSE
    if( n < 5 )
        LogicError("Must be at least 5x5 to have a fourth-order symbol");
    typedef Complex<Real> C;
    Zeros( A, n, n );
    FillDiagonal( A, C(-1),     4 );
    FillDiagonal( A, C(-3,-2),  3 );
    FillDiagonal( A, C( 0, 1),  2 );
    FillDiagonal( A, C(1),      1 );
    FillDiagonal( A, C(10),    -1 );
    FillDiagonal( A, C( 3, 1), -2 );
    FillDiagonal( A, C(4),     -3 );
    FillDiagonal( A, C( 0, 1), -4 );
}

#define PROTO(Real) \
  template void Whale( Matrix<Complex<Real>>& A, Int n ); \
  template void Whale( AbstractDistMatrix<Complex<Real>>& A, Int n );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_QUAD
#include <El/macros/Instantiate.h>

} // namespace El
