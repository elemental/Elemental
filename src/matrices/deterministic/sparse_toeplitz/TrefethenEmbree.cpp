/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>
#include <El/matrices.hpp>

namespace El {

// The symbol of the matrix described at the beginning of Chapter II of 
// Lloyd N. Trefethen and Mark Embree's "Spectra and Pseudospectra" is
//     f(z) = 2 z^{-3} - z^{-2} + 2i z^{-1} - 4 z^2 - 2i z^3.
// For lack of a better name, we will refer to such a Toeplitz matrix as a 
// Trefethen-Embree matrix.

template<typename Real> 
void TrefethenEmbree( Matrix<Complex<Real>>& A, Int n )
{
    EL_DEBUG_CSE
    if( n < 4 )
        LogicError("Must be at least 4x4 to have a third-order symbol");
    typedef Complex<Real> C;
    Zeros( A, n, n );
    FillDiagonal( A, C(2),     3 );
    FillDiagonal( A, C(-1),    2 );
    FillDiagonal( A, C(0,2),   1 );
    FillDiagonal( A, C(-4),   -2 );
    FillDiagonal( A, C(0,-2), -3 );
}

template<typename Real>
void TrefethenEmbree( AbstractDistMatrix<Complex<Real>>& A, Int n )
{
    EL_DEBUG_CSE
    if( n < 4 )
        LogicError("Must be at least 4x4 to have a third-order symbol");
    typedef Complex<Real> C;
    Zeros( A, n, n );
    FillDiagonal( A, C(2),     3 );
    FillDiagonal( A, C(-1),    2 );
    FillDiagonal( A, C(0,2),   1 );
    FillDiagonal( A, C(-4),   -2 );
    FillDiagonal( A, C(0,-2), -3 );
}

#define PROTO(Real) \
  template void TrefethenEmbree( Matrix<Complex<Real>>& A, Int n ); \
  template void TrefethenEmbree \
  ( AbstractDistMatrix<Complex<Real>>& A, Int n );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_QUAD
#include <El/macros/Instantiate.h>

} // namespace El
