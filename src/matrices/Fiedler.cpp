/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F> 
void Fiedler( Matrix<F>& A, const std::vector<F>& c )
{
    DEBUG_ONLY(CallStackEntry cse("Fiedler"))
    const Int n = c.size();
    A.Resize( n, n );
    IndexDependentFill( A, [&]( Int i, Int j ) { return Abs(c[i]-c[j]); } );
}

template<typename F>
void Fiedler( AbstractDistMatrix<F>& A, const std::vector<F>& c )
{
    DEBUG_ONLY(CallStackEntry cse("Fiedler"))
    const Int n = c.size();
    A.Resize( n, n );
    IndexDependentFill( A, [&]( Int i, Int j ) { return Abs(c[i]-c[j]); } );
}

template<typename F>
void Fiedler( AbstractBlockDistMatrix<F>& A, const std::vector<F>& c )
{
    DEBUG_ONLY(CallStackEntry cse("Fiedler"))
    const Int n = c.size();
    A.Resize( n, n );
    IndexDependentFill( A, [&]( Int i, Int j ) { return Abs(c[i]-c[j]); } );
}

#define PROTO(F) \
  template void Fiedler( Matrix<F>& A, const std::vector<F>& c ); \
  template void Fiedler( AbstractDistMatrix<F>& A, const std::vector<F>& c ); \
  template void Fiedler \
  ( AbstractBlockDistMatrix<F>& A, const std::vector<F>& c );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
