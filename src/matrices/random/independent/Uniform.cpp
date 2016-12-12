/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>
#include <El/matrices.hpp>

namespace El {

// Draw each entry from a uniform PDF over a closed ball.

template<typename T>
void MakeUniform( Matrix<T>& A, T center, Base<T> radius )
{
    EL_DEBUG_CSE
    auto sampleBall = [=]() { return SampleBall(center,radius); };
    EntrywiseFill( A, function<T()>(sampleBall) );
}

template<typename T>
void Uniform( Matrix<T>& A, Int m, Int n, T center, Base<T> radius )
{
    EL_DEBUG_CSE
    A.Resize( m, n );
    MakeUniform( A, center, radius );
}

template<typename T>
void MakeUniform( AbstractDistMatrix<T>& A, T center, Base<T> radius )
{
    EL_DEBUG_CSE
    if( A.RedundantRank() == 0 )
        MakeUniform( A.Matrix(), center, radius );
    Broadcast( A, A.RedundantComm(), 0 );
}

template<typename T>
void Uniform( AbstractDistMatrix<T>& A, Int m, Int n, T center, Base<T> radius )
{
    EL_DEBUG_CSE
    A.Resize( m, n );
    MakeUniform( A, center, radius );
}

template<typename T>
void MakeUniform( DistMultiVec<T>& X, T center, Base<T> radius )
{
    EL_DEBUG_CSE
    const int localHeight = X.LocalHeight();
    const int width = X.Width();
    for( int j=0; j<width; ++j )
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
            X.SetLocal( iLocal, j, SampleBall(center,radius) );
}

template<typename T>
void Uniform( DistMultiVec<T>& A, Int m, Int n, T center, Base<T> radius )
{
    EL_DEBUG_CSE
    A.Resize( m, n );
    MakeUniform( A, center, radius );
}

#define PROTO(T) \
  template void MakeUniform \
  ( Matrix<T>& A, T center, Base<T> radius ); \
  template void MakeUniform \
  ( AbstractDistMatrix<T>& A, T center, Base<T> radius ); \
  template void MakeUniform( DistMultiVec<T>& A, T center, Base<T> radius ); \
  template void Uniform \
  ( Matrix<T>& A, Int m, Int n, T center, Base<T> radius ); \
  template void Uniform \
  ( AbstractDistMatrix<T>& A, Int m, Int n, T center, Base<T> radius ); \
  template void Uniform \
  ( DistMultiVec<T>& A, Int m, Int n, T center, Base<T> radius );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
