/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T> 
void Hankel( Matrix<T>& A, Int m, Int n, const std::vector<T>& a )
{
    DEBUG_ONLY(CallStackEntry cse("Hankel"))
    const Int length = m+n-1;
    if( a.size() != (Unsigned)length )
        LogicError("a was the wrong size");
    A.Resize( m, n );
    IndexDependentFill( A, [&]( Int i, Int j ) { return a[i+j]; } );
}

template<typename T>
void Hankel( AbstractDistMatrix<T>& A, Int m, Int n, const std::vector<T>& a )
{
    DEBUG_ONLY(CallStackEntry cse("Hankel"))
    const Int length = m+n-1;
    if( a.size() != (Unsigned)length )
        LogicError("a was the wrong size");
    A.Resize( m, n );
    IndexDependentFill( A, [&]( Int i, Int j ) { return a[i+j]; } );
}

template<typename T>
void Hankel
( AbstractBlockDistMatrix<T>& A, Int m, Int n, const std::vector<T>& a )
{
    DEBUG_ONLY(CallStackEntry cse("Hankel"))
    const Int length = m+n-1;
    if( a.size() != (Unsigned)length )
        LogicError("a was the wrong size");
    A.Resize( m, n );
    IndexDependentFill( A, [&]( Int i, Int j ) { return a[i+j]; } );
}

#define PROTO(T) \
  template void Hankel \
  ( Matrix<T>& A, Int m, Int n, const std::vector<T>& a ); \
  template void Hankel \
  ( AbstractDistMatrix<T>& A, Int m, Int n, const std::vector<T>& a ); \
  template void Hankel \
  ( AbstractBlockDistMatrix<T>& A, Int m, Int n, const std::vector<T>& a );

PROTO(Int)
PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
