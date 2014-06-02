/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include "El/matrices/Uniform.hpp"

namespace El {

// Draw each entry from a uniform PDF over a closed ball.

template<typename T>
void MakeUniform( Matrix<T>& A, T center, Base<T> radius )
{
    DEBUG_ONLY(CallStackEntry cse("MakeUniform"))
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, SampleBall( center, radius ) );
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
    if( A.RedundantSize() == 1 )
    {
        MakeUniform( A.Matrix(), center, radius );
    }
    else if( A.Participating() && A.LocalHeight() == A.LDim() )
    {
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        if( A.RedundantRank() == 0 )
            MakeUniform( A.Matrix(), center, radius );
        mpi::Broadcast
        ( A.Buffer(), localHeight*localWidth, 0, A.RedundantComm() );
    }
    else if( A.Participating() )
    {
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        const Int bufSize = localHeight*localWidth;   
        std::vector<T> buffer( bufSize );
        if( A.RedundantRank() == 0 )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    buffer[iLoc+jLoc*localHeight] = 
                        SampleBall( center, radius );
        }
        mpi::Broadcast( buffer.data(), bufSize, 0, A.RedundantComm() );
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( A.Buffer(0,jLoc), &buffer[jLoc*localHeight], localHeight );
    }
}

template<typename T>
void MakeUniform( AbstractBlockDistMatrix<T>& A, T center, Base<T> radius )
{
    DEBUG_ONLY(CallStackEntry cse("MakeUniform"))
    if( A.RedundantSize() == 1 )
    {
        MakeUniform( A.Matrix(), center, radius );
    }
    else if( A.Participating() && A.LocalHeight() == A.LDim() )
    {
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        if( A.RedundantRank() == 0 )
            MakeUniform( A.Matrix(), center, radius );
        mpi::Broadcast
        ( A.Buffer(), localHeight*localWidth, 0, A.RedundantComm() );
    }
    else if( A.Participating() )
    {
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        const Int bufSize = localHeight*localWidth;   
        std::vector<T> buffer( bufSize );
        if( A.RedundantRank() == 0 )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                    buffer[iLoc+jLoc*localHeight] = 
                        SampleBall( center, radius );
        }
        mpi::Broadcast( buffer.data(), bufSize, 0, A.RedundantComm() );
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( A.Buffer(0,jLoc), &buffer[jLoc*localHeight], localHeight );
    }
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
