/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename S,typename T> 
void Toeplitz( Matrix<S>& A, Int m, Int n, const vector<T>& a )
{
    DEBUG_ONLY(CallStackEntry cse("Toeplitz"))
    const Int length = m+n-1;
    if( a.size() != Unsigned(length) )
        LogicError("a was the wrong size");
    A.Resize( m, n );
    auto toeplitzFill = [&]( Int i, Int j ) { return a[i-j+(n-1)]; };
    IndexDependentFill( A, function<S(Int,Int)>(toeplitzFill) );
}

template<typename S,typename T>
void Toeplitz( AbstractDistMatrix<S>& A, Int m, Int n, const vector<T>& a )
{
    DEBUG_ONLY(CallStackEntry cse("Toeplitz"))
    const Int length = m+n-1;
    if( a.size() != Unsigned(length) )
        LogicError("a was the wrong size");
    A.Resize( m, n );
    auto toeplitzFill = [&]( Int i, Int j ) { return a[i-j+(n-1)]; };
    IndexDependentFill( A, function<S(Int,Int)>(toeplitzFill) );
}

template<typename S,typename T>
void Toeplitz
( AbstractBlockDistMatrix<S>& A, Int m, Int n, const vector<T>& a )
{
    DEBUG_ONLY(CallStackEntry cse("Toeplitz"))
    const Int length = m+n-1;
    if( a.size() != Unsigned(length) )
        LogicError("a was the wrong size");
    A.Resize( m, n );
    auto toeplitzFill = [&]( Int i, Int j ) { return a[i-j+(n-1)]; };
    IndexDependentFill( A, function<S(Int,Int)>(toeplitzFill) );
}

#define PROTO_TYPES(T1,T2) \
  template void Toeplitz \
  ( Matrix<T1>& A, \
    const Int m, const Int n, const vector<T2>& a ); \
  template void Toeplitz \
  ( AbstractDistMatrix<T1>& A, \
    const Int m, const Int n, const vector<T2>& a ); \
  template void Toeplitz \
  ( AbstractBlockDistMatrix<T1>& A, \
    const Int m, const Int n, const vector<T2>& a );

#define PROTO_INT(T) PROTO_TYPES(T,T)

#define PROTO_REAL(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,T)

#define PROTO_COMPLEX(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,Base<T>) \
  PROTO_TYPES(T,T)

#include "El/macros/Instantiate.h"

} // namespace El
