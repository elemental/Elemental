/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_UNIFORM_HPP
#define ELEM_UNIFORM_HPP

namespace elem {

// Draw each entry from a uniform PDF over a closed ball.
template<typename T>
inline void
MakeUniform( Matrix<T>& A, T center=0, Base<T> radius=1 )
{
    DEBUG_ONLY(CallStackEntry cse("MakeUniform"))
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, SampleBall( center, radius ) );
}

template<typename T>
inline void
Uniform( Matrix<T>& A, Int m, Int n, T center=0, Base<T> radius=1 )
{
    DEBUG_ONLY(CallStackEntry cse("Uniform"))
    A.Resize( m, n );
    MakeUniform( A, center, radius );
}

template<typename T,Dist U,Dist V>
inline void
MakeUniform( DistMatrix<T,U,V>& A, T center=0, Base<T> radius=1 )
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

template<typename T,Dist U,Dist V>
inline void
MakeUniform( BlockDistMatrix<T,U,V>& A, T center=0, Base<T> radius=1 )
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

template<typename T,Dist U,Dist V>
inline void
Uniform( DistMatrix<T,U,V>& A, Int m, Int n, T center=0, Base<T> radius=1 )
{
    DEBUG_ONLY(CallStackEntry cse("Uniform"))
    A.Resize( m, n );
    MakeUniform( A, center, radius );
}

template<typename T,Dist U,Dist V>
inline void
Uniform( BlockDistMatrix<T,U,V>& A, Int m, Int n, T center=0, Base<T> radius=1 )
{
    DEBUG_ONLY(CallStackEntry cse("Uniform"))
    A.Resize( m, n );
    MakeUniform( A, center, radius );
}

} // namespace elem

#endif // ifndef ELEM_UNIFORM_HPP
