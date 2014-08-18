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
void Lauchli( Matrix<T>& A, Int n, T mu )
{
    DEBUG_ONLY(CallStackEntry cse("Lauchli"))
    A.Resize( n+1, n );

    auto ABlock = A( IR(0,1), IR(0,n) );
    Fill( ABlock, T(1) );

    std::vector<T> d(n,mu);
    ABlock = A( IR(1,n+1), IR(0,n) );
    Diagonal( ABlock, d );
}

template<typename T>
void Lauchli( AbstractDistMatrix<T>& A, Int n, T mu )
{
    DEBUG_ONLY(CallStackEntry cse("Lauchli"))
    A.Resize( n+1, n );

    std::unique_ptr<AbstractDistMatrix<T>> ABlock( A.Construct(A.Grid()) );
    View( *ABlock, A, IR(0,1), IR(0,n) );
    Fill( *ABlock, T(1) );

    std::vector<T> d(n,mu);
    View( *ABlock, A, IR(1,n+1), IR(0,n) );
    Diagonal( *ABlock, d );
}

#define PROTO(T) \
  template void Lauchli( Matrix<T>& A, Int n, T mu ); \
  template void Lauchli( AbstractDistMatrix<T>& A, Int n, T mu );

#include "El/macros/Instantiate.h"

} // namespace El
