/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void RealPart( const Matrix<T>& A, Matrix<Base<T>>& AReal )
{
    DEBUG_ONLY(CallStackEntry cse("EntrywiseMap"))
    const Int m = A.Height();
    const Int n = A.Width();
    AReal.Resize( m, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            AReal.Set( i, j, A.GetRealPart(i,j) );
}

template<typename T>
void RealPart
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<Base<T>>& AReal )
{ 
    auto realPart = []( T alpha ) { return RealPart(alpha); };
    std::function<Base<T>(T)> realLambda( realPart );
    EntrywiseMap( A, AReal, realLambda );
}

template<typename T>
void RealPart
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<Base<T>>& AReal )
{ 
    auto realPart = []( T alpha ) { return RealPart(alpha); };
    std::function<Base<T>(T)> realLambda( realPart );
    EntrywiseMap( A, AReal, realLambda );
}

#define PROTO(T) \
  template void RealPart( const Matrix<T>& A, Matrix<Base<T>>& AReal ); \
  template void RealPart \
  ( const AbstractDistMatrix<T>& A, AbstractDistMatrix<Base<T>>& AReal ); \
  template void RealPart \
  ( const AbstractBlockDistMatrix<T>& A, \
          AbstractBlockDistMatrix<Base<T>>& AReal );

#include "El/macros/Instantiate.h"

} // namespace El
