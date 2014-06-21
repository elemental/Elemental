/*
   Copyright (c) 2009-2014, Jack Poulson
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
    EntrywiseFill( A, [=]() { return SampleBall( center, radius ); } );
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
    A.BroadcastOver( A.RedundantComm(), 0 );
}

template<typename T>
void MakeUniform( AbstractBlockDistMatrix<T>& A, T center, Base<T> radius )
{
    DEBUG_ONLY(CallStackEntry cse("MakeUniform"))
    if( A.RedundantRank() == 0 )
        MakeUniform( A.Matrix(), center, radius );
    A.BroadcastOver( A.RedundantComm(), 0 );
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

#define PROTO(T) \
  template void MakeUniform \
  ( Matrix<T>& A, T center, Base<T> radius ); \
  template void MakeUniform \
  ( AbstractDistMatrix<T>& A, T center, Base<T> radius ); \
  template void MakeUniform \
  ( AbstractBlockDistMatrix<T>& A, T center, Base<T> radius ); \
  template void Uniform \
  ( Matrix<T>& A, Int m, Int n, T center, Base<T> radius ); \
  template void Uniform \
  ( AbstractDistMatrix<T>& A, Int m, Int n, T center, Base<T> radius ); \
  template void Uniform \
  ( AbstractBlockDistMatrix<T>& A, Int m, Int n, T center, Base<T> radius )

PROTO(Int);
#ifndef EL_DISABLE_FLOAT
PROTO(float);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<float>);
#endif // ifndef EL_DISABLE_COMPLEX
#endif // ifndef EL_DISABLE_FLOAT

PROTO(double);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<double>);
#endif // ifndef EL_DISABLE_COMPLEX

} // namespace El
