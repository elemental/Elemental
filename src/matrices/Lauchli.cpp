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
void Lauchli( Matrix<T>& A, Int n, T mu )
{
    DEBUG_ONLY(CallStackEntry cse("Lauchli"))
    Zeros( A, n+1, n );

    // Set the first row to all ones
    auto a0 = A( IR(0,1), IR(0,n) );
    Fill( a0, T(1) );

    // Set the subdiagonal to mu
    FillDiagonal( A, mu, -1 );
}

template<typename T>
void Lauchli( AbstractDistMatrix<T>& A, Int n, T mu )
{
    DEBUG_ONLY(CallStackEntry cse("Lauchli"))
    Zeros( A, n+1, n );

    // Set the first row to all ones
    unique_ptr<AbstractDistMatrix<T>> a0( A.Construct(A.Grid(),A.Root()) );
    View( *a0, A, IR(0,1), IR(0,n) );
    Fill( *a0, T(1) );

    // Set the subdiagonal to mu
    FillDiagonal( A, mu, -1 );
}

#define PROTO(T) \
  template void Lauchli( Matrix<T>& A, Int n, T mu ); \
  template void Lauchli( AbstractDistMatrix<T>& A, Int n, T mu );

#include "El/macros/Instantiate.h"

} // namespace El
