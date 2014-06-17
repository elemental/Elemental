/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename S,typename T> 
void Diagonal( Matrix<S>& D, const std::vector<T>& d )
{
    DEBUG_ONLY(CallStackEntry cse("Diagonal"))
    const Int n = d.size();
    Zeros( D, n, n );

    for( Int j=0; j<n; ++j )
        D.Set( j, j, d[j] );
}

template<typename S,typename T>
void Diagonal( AbstractDistMatrix<S>& D, const std::vector<T>& d )
{
    DEBUG_ONLY(CallStackEntry cse("Diagonal"))
    const Int n = d.size();
    Zeros( D, n, n );

    const Int localWidth = D.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = D.GlobalCol(jLoc);
        D.Set( j, j, d[j] );
    }
}

template<typename S,typename T>
void Diagonal( AbstractBlockDistMatrix<S>& D, const std::vector<T>& d )
{
    DEBUG_ONLY(CallStackEntry cse("Diagonal"))
    const Int n = d.size();
    Zeros( D, n, n );

    const Int localWidth = D.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = D.GlobalCol(jLoc);
        D.Set( j, j, d[j] );
    }
}

#define PROTO_TYPES(S,T) \
  template void Diagonal( Matrix<S>& D, const std::vector<T>& d ); \
  template void Diagonal( AbstractDistMatrix<S>& D, const std::vector<T>& d ); \
  template void Diagonal \
  ( AbstractBlockDistMatrix<S>& D, const std::vector<T>& d );

#define PROTO_INT(S) PROTO_TYPES(S,S)

#define PROTO_REAL(S) \
  PROTO_TYPES(S,Int) \
  PROTO_TYPES(S,S)

#define PROTO_CPX(S) \
  PROTO_TYPES(S,Int) \
  PROTO_TYPES(S,Base<S>) \
  PROTO_TYPES(S,S)

PROTO_INT(Int)
PROTO_REAL(float)
PROTO_REAL(double)
PROTO_CPX(Complex<float>)
PROTO_CPX(Complex<double>)

} // namespace El
