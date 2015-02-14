/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Draw each entry from a uniform PDF over a closed ball.

template<typename T>
void MakeUniform( Matrix<T>& A, T center, Base<T> radius )
{
    DEBUG_ONLY(CallStackEntry cse("MakeUniform"))
    auto sampleBall = [=]() { return SampleBall(center,radius); };
    EntrywiseFill( A, function<T()>(sampleBall) );
}

template<typename T>
void Uniform( Matrix<T>& A, Int m, Int n, T center, Base<T> radius )
{
    DEBUG_ONLY(CallStackEntry cse("Uniform"))
    A.Resize( m, n );
    MakeUniform( A, center, radius );
}

template<typename T>
void MakeUniform( AbstractDistMatrix<T>& A, T center, Base<T> radius )
{
    DEBUG_ONLY(CallStackEntry cse("MakeUniform"))
    if( A.RedundantRank() == 0 )
        MakeUniform( A.Matrix(), center, radius );
    Broadcast( A, A.RedundantComm(), 0 );
}

template<typename T>
void MakeUniform( AbstractBlockDistMatrix<T>& A, T center, Base<T> radius )
{
    DEBUG_ONLY(CallStackEntry cse("MakeUniform"))
    if( A.RedundantRank() == 0 )
        MakeUniform( A.Matrix(), center, radius );
    Broadcast( A, A.RedundantComm(), 0 );
}

template<typename T>
void Uniform( AbstractDistMatrix<T>& A, Int m, Int n, T center, Base<T> radius )
{
    DEBUG_ONLY(CallStackEntry cse("Uniform"))
    A.Resize( m, n );
    MakeUniform( A, center, radius );
}

template<typename T>
void Uniform
( AbstractBlockDistMatrix<T>& A, Int m, Int n, T center, Base<T> radius )
{
    DEBUG_ONLY(CallStackEntry cse("Uniform"))
    A.Resize( m, n );
    MakeUniform( A, center, radius );
}

template<typename T>
void MakeUniform( DistMultiVec<T>& X, T center, Base<T> radius )
{
    DEBUG_ONLY(CallStackEntry cse("MakeUniform"))
    const int localHeight = X.LocalHeight();
    const int width = X.Width();
    for( int j=0; j<width; ++j )
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            X.SetLocal( iLocal, j, SampleBall(center,radius) );
}

template<typename T>
void Uniform( DistMultiVec<T>& A, Int m, Int n, T center, Base<T> radius )
{
    DEBUG_ONLY(CallStackEntry cse("Uniform"))
    A.Resize( m, n );
    MakeUniform( A, center, radius );
}

#define PROTO(T) \
  template void MakeUniform \
  ( Matrix<T>& A, T center, Base<T> radius ); \
  template void MakeUniform \
  ( AbstractDistMatrix<T>& A, T center, Base<T> radius ); \
  template void MakeUniform \
  ( AbstractBlockDistMatrix<T>& A, T center, Base<T> radius ); \
  template void MakeUniform( DistMultiVec<T>& A, T center, Base<T> radius ); \
  template void Uniform \
  ( Matrix<T>& A, Int m, Int n, T center, Base<T> radius ); \
  template void Uniform \
  ( AbstractDistMatrix<T>& A, Int m, Int n, T center, Base<T> radius ); \
  template void Uniform \
  ( AbstractBlockDistMatrix<T>& A, Int m, Int n, T center, Base<T> radius ); \
  template void Uniform \
  ( DistMultiVec<T>& A, Int m, Int n, T center, Base<T> radius );

#include "El/macros/Instantiate.h"

} // namespace El
