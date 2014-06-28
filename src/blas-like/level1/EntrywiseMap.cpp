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
void EntrywiseMap( Matrix<T>& A, std::function<T(T)> func )
{
    DEBUG_ONLY(CallStackEntry cse("EntrywiseMap"))
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, func(A.Get(i,j)) );
}

template<typename T>
void EntrywiseMap( AbstractDistMatrix<T>& A, std::function<T(T)> func )
{ EntrywiseMap( A.Matrix(), func ); }

template<typename T>
void EntrywiseMap( AbstractBlockDistMatrix<T>& A, std::function<T(T)> func )
{ EntrywiseMap( A.Matrix(), func ); }

#define PROTO(T) \
  template void EntrywiseMap( Matrix<T>& A, std::function<T(T)> func ); \
  template void EntrywiseMap \
  ( AbstractDistMatrix<T>& A, std::function<T(T)> func ); \
  template void EntrywiseMap \
  ( AbstractBlockDistMatrix<T>& A, std::function<T(T)> func );

PROTO(Int)
PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
