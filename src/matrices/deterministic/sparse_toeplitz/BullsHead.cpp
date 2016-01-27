/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// The "bull's head matrix" is defined to have the symbol:
//   f(z) = 2i z^{-1} + z^2 + 7/10 z^3.
// Please see 
//   L. Reichel and L. N. Trefethen, "Eigenvalues and pseudo-eigenvalues of 
//   Toeplitz matrices", Linear Algebra Appl., 1992.

template<typename Real> 
void BullsHead( Matrix<Complex<Real>>& A, Int n )
{
    DEBUG_ONLY(CSE cse("BullsHead"))
    if( n < 4 )
        LogicError("Must be at least 4x4 to have a third-order symbol");
    typedef Complex<Real> C;
    Zeros( A, n, n );
    FillDiagonal( A, C(0,2),      1 );
    FillDiagonal( A, C(1),       -2 );
    FillDiagonal( A, C(7)/C(10), -3 );
}

template<typename Real>
void BullsHead( AbstractDistMatrix<Complex<Real>>& A, Int n )
{
    DEBUG_ONLY(CSE cse("BullsHead"))
    if( n < 4 )
        LogicError("Must be at least 4x4 to have a third-order symbol");
    typedef Complex<Real> C;
    Zeros( A, n, n );
    FillDiagonal( A, C(0,2),      1 );
    FillDiagonal( A, C(1),       -2 );
    FillDiagonal( A, C(7)/C(10), -3 );
}

#define PROTO(Real) \
   template void BullsHead( Matrix<Complex<Real>>& A, Int n ); \
   template void BullsHead( AbstractDistMatrix<Complex<Real>>& A, Int n );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
