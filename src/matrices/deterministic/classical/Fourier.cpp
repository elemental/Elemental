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
void Fourier( Matrix<Complex<Real>>& A, Int n )
{
    DEBUG_CSE
    A.Resize( n, n );
    const Real pi = 4*Atan( Real(1) );
    const Real nSqrt = Sqrt( Real(n) );
    auto fourierFill = 
      [=]( Int i, Int j ) -> Complex<Real>
      { const Real theta = -2*pi*i*j/n;
        return Complex<Real>(Cos(theta),Sin(theta))/nSqrt; };
    IndexDependentFill( A, function<Complex<Real>(Int,Int)>(fourierFill) );
}

template<typename Real>
void Fourier( AbstractDistMatrix<Complex<Real>>& A, Int n )
{
    DEBUG_CSE
    A.Resize( n, n );
    const Real pi = 4*Atan( Real(1) );
    const Real nSqrt = Sqrt( Real(n) );
    auto fourierFill = 
      [=]( Int i, Int j ) -> Complex<Real>
      { const Real theta = -2*pi*i*j/n;
        return Complex<Real>(Cos(theta),Sin(theta))/nSqrt; };
    IndexDependentFill( A, function<Complex<Real>(Int,Int)>(fourierFill) );
}

#define PROTO(Real) \
  template void Fourier( Matrix<Complex<Real>>& A, Int n ); \
  template void Fourier( AbstractDistMatrix<Complex<Real>>& A, Int n );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_QUAD
#include <El/macros/Instantiate.h>

} // namespace El
