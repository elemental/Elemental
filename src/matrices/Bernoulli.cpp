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
void Bernoulli( Matrix<T>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Bernoulli"))
    A.Resize( m, n );
    auto coinFlip = [=]() { return T(CoinFlip()); };
    EntrywiseFill( A, std::function<T()>(coinFlip) );
}

template<typename T>
void Bernoulli( AbstractDistMatrix<T>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Bernoulli"))
    A.Resize( m, n );
    if( A.RedundantRank() == 0 )
        Bernoulli( A.Matrix(), A.LocalHeight(), A.LocalWidth() );
    A.BroadcastOver( A.RedundantComm(), 0 );
}

template<typename T>
void Bernoulli( AbstractBlockDistMatrix<T>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Bernoulli"))
    A.Resize( m, n );
    if( A.RedundantRank() == 0 )
        Bernoulli( A.Matrix(), A.LocalHeight(), A.LocalWidth() );
    A.BroadcastOver( A.RedundantComm(), 0 );
}


#define PROTO(T) \
  template void Bernoulli( Matrix<T>& A, Int m, Int n ); \
  template void Bernoulli( AbstractDistMatrix<T>& A, Int m, Int n ); \
  template void Bernoulli( AbstractBlockDistMatrix<T>& A, Int m, Int n );

#include "El/macros/Instantiate.h"

} // namespace El
