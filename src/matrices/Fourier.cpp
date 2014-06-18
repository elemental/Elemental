/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename Real>
void Fourier( Matrix<Complex<Real>>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Fourier"))
    A.Resize( n, n );
    const Real pi = 4*Atan( Real(1) );
    const Real nSqrt = Sqrt( Real(n) );
    IndexDependentFill
    ( A, [=]( Int i, Int j ) 
         { 
             const Real theta = -2*pi*i*j/n;
             return Complex<Real>(Cos(theta),Sin(theta))/nSqrt;
         } );
}

template<typename Real>
void Fourier( AbstractDistMatrix<Complex<Real>>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Fourier"))
    A.Resize( n, n );
    const Real pi = 4*Atan( Real(1) );
    const Real nSqrt = Sqrt( Real(n) );
    IndexDependentFill
    ( A, [=]( Int i, Int j ) 
         { 
             const Real theta = -2*pi*i*j/n;
             return Complex<Real>(Cos(theta),Sin(theta))/nSqrt;
         } );
}

template<typename Real>
void Fourier( AbstractBlockDistMatrix<Complex<Real>>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Fourier"))
    A.Resize( n, n );
    const Real pi = 4*Atan( Real(1) );
    const Real nSqrt = Sqrt( Real(n) );
    IndexDependentFill
    ( A, [=]( Int i, Int j ) 
         { 
             const Real theta = -2*pi*i*j/n;
             return Complex<Real>(Cos(theta),Sin(theta))/nSqrt;
         } );
}

#define PROTO(Real) \
  template void Fourier( Matrix<Complex<Real>>& A, Int n ); \
  template void Fourier( AbstractDistMatrix<Complex<Real>>& A, Int n ); \
  template void Fourier( AbstractBlockDistMatrix<Complex<Real>>& A, Int n );

PROTO(float)
PROTO(double)

} // namespace El
