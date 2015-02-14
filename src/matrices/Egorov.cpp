/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename Real>
void Egorov
( Matrix<Complex<Real>>& A, function<Real(Int,Int)> phase, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Egorov"))
    A.Resize( n, n );
    auto egorovFill = 
      [&]( Int i, Int j ) -> Complex<Real>
      { const Real theta = phase(i,j);
        return Complex<Real>(Cos(theta),Sin(theta)); }; 
    IndexDependentFill( A, function<Complex<Real>(Int,Int)>(egorovFill) );
}

template<typename Real>
void Egorov
( AbstractDistMatrix<Complex<Real>>& A, 
  function<Real(Int,Int)> phase, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Egorov"))
    A.Resize( n, n );
    auto egorovFill = 
      [&]( Int i, Int j ) -> Complex<Real>
      { const Real theta = phase(i,j);
        return Complex<Real>(Cos(theta),Sin(theta)); }; 
    IndexDependentFill( A, function<Complex<Real>(Int,Int)>(egorovFill) );
}

template<typename Real>
void Egorov
( AbstractBlockDistMatrix<Complex<Real>>& A,
  function<Real(Int,Int)> phase, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Egorov"))
    A.Resize( n, n );
    auto egorovFill = 
      [&]( Int i, Int j ) -> Complex<Real>
      { const Real theta = phase(i,j);
        return Complex<Real>(Cos(theta),Sin(theta)); }; 
    IndexDependentFill( A, function<Complex<Real>(Int,Int)>(egorovFill) );
}

#define PROTO(Real) \
  template void Egorov \
  ( Matrix<Complex<Real>>& A, function<Real(Int,Int)> phase, Int n ); \
  template void Egorov \
  ( AbstractDistMatrix<Complex<Real>>& A, \
    function<Real(Int,Int)> phase, Int n ); \
  template void Egorov \
  ( AbstractBlockDistMatrix<Complex<Real>>& A, \
    function<Real(Int,Int)> phase, Int n );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
