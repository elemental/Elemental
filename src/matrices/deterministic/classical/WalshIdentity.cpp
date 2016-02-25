/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void WalshIdentity( Matrix<T>& A, Int k, bool binary )
{
    DEBUG_ONLY(CSE cse("WalshIdentity"))
    if( k < 1 )
      LogicError("Walsh matrices are only defined for k>=1");
    const Unsigned n = 1u<<k;
    A.Resize( n, 2*n );
    auto AL = A( IR(0,n), IR(0,n) );
    auto AR = A( IR(0,n), IR(n,2*n) );
    Walsh( AL, k, binary );
    Identity( AR, n, n );
}

template<typename T>
void WalshIdentity( ElementalMatrix<T>& A, Int k, bool binary )
{
    DEBUG_ONLY(CSE cse("WalshIdentity"))
    if( k < 1 )
      LogicError("Walsh matrices are only defined for k>=1");
    const Unsigned n = 1u<<k;
    A.Resize( n, 2*n );
    unique_ptr<ElementalMatrix<T>> AL( A.Construct(A.Grid(),A.Root()) );
    unique_ptr<ElementalMatrix<T>> AR( A.Construct(A.Grid(),A.Root()) );
    View( *AL, A, IR(0,n), IR(0,n) );
    View( *AR, A, IR(0,n), IR(n,2*n) );
    Walsh( *AL, k, binary );
    Identity( *AR, n, n );
}

#define PROTO(T) \
  template void WalshIdentity \
  ( Matrix<T>& A, Int k, bool binary ); \
  template void WalshIdentity( ElementalMatrix<T>& A, Int k, bool binary );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
