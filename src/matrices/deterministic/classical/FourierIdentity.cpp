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

template<typename Real>
void FourierIdentity( Matrix<Complex<Real>>& A, Int n )
{
    DEBUG_CSE
    A.Resize( n, 2*n );
    auto AL = A( IR(0,n), IR(0,n) );
    auto AR = A( IR(0,n), IR(n,2*n) );
    Fourier( AL, n );
    Identity( AR, n, n );
}

template<typename Real>
void FourierIdentity( ElementalMatrix<Complex<Real>>& A, Int n )
{
    DEBUG_CSE
    typedef Complex<Real> C;

    A.Resize( n, 2*n );
    unique_ptr<ElementalMatrix<C>> AL( A.Construct(A.Grid(),A.Root()) );
    unique_ptr<ElementalMatrix<C>> AR( A.Construct(A.Grid(),A.Root()) );
    View( *AL, A, IR(0,n), IR(0,n) );
    View( *AR, A, IR(0,n), IR(n,2*n) );
    Fourier( *AL, n );
    Identity( *AR, n, n );
}

#define PROTO(Real) \
  template void FourierIdentity \
  ( Matrix<Complex<Real>>& A, Int n ); \
  template void FourierIdentity \
  ( ElementalMatrix<Complex<Real>>& A, Int n );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_QUAD
#include <El/macros/Instantiate.h>

} // namespace El
