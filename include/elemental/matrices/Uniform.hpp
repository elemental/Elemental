/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_UNIFORM_HPP
#define ELEM_MATRICES_UNIFORM_HPP

namespace elem {

// Draw each entry from a uniform PDF over a closed ball.
template<typename T>
inline void
MakeUniform( Matrix<T>& A, T center=0, BASE(T) radius=1 )
{
#ifndef RELEASE
    CallStackEntry cse("MakeUniform");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, SampleBall( center, radius ) );
}

template<typename T>
inline void
Uniform( Matrix<T>& A, Int m, Int n, T center=0, BASE(T) radius=1 )
{
#ifndef RELEASE
    CallStackEntry cse("Uniform");
#endif
    A.ResizeTo( m, n );
    MakeUniform( A, center, radius );
}

template<typename T>
inline Matrix<T>
Uniform( Int m, Int n, T center=0, BASE(T) radius=1 )
{
    Matrix<T> A( m, n );
    MakeUniform( A, center, radius );
    return A;
}

namespace internal {

template<typename T,Distribution U,Distribution V>
struct MakeUniformHelper
{
    static void Func( DistMatrix<T,U,V>& A, T center, BASE(T) radius );  
};

template<typename T>
struct MakeUniformHelper<T,CIRC,CIRC>
{
    static void Func( DistMatrix<T,CIRC,CIRC>& A, T center, BASE(T) radius )
    {
        if( A.Grid().VCRank() == A.Root() )
        {
            const Int height = A.Height(); 
            const Int width = A.Width();
            for( Int j=0; j<width; ++j )
                for( Int i=0; i<height; ++i )
                    A.SetLocal( i, j, SampleBall( center, radius ) );
        }
    }
};

template<typename T>
struct MakeUniformHelper<T,MC,MR>
{
    static void Func( DistMatrix<T,MC,MR>& A, T center, BASE(T) radius )
    {
        const Int localHeight = A.LocalHeight(); 
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                A.SetLocal( iLoc, jLoc, SampleBall( center, radius ) );
    }
};

template<typename T>
struct MakeUniformHelper<T,MC,STAR>
{
    static void Func( DistMatrix<T,MC,STAR>& A, T center, BASE(T) radius )
    {
        const Grid& grid = A.Grid();
        if( grid.InGrid() )
        {
            const Int n = A.Width();
            const Int localHeight = A.LocalHeight();
            const Int bufSize = localHeight*n;
            std::vector<T> buffer( bufSize );

            // Create random matrix on process column 0, then broadcast
            if( grid.Col() == 0 )
            {
                for( Int j=0; j<n; ++j )
                    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                        buffer[iLoc+j*localHeight] = 
                            SampleBall( center, radius );
            }
            mpi::Broadcast( buffer.data(), bufSize, 0, grid.RowComm() );

            // Unpack
            T* localBuffer = A.Buffer();
            const Int ldim = A.LDim();
            PARALLEL_FOR
            for( Int j=0; j<n; ++j )
            {
                const T* bufferCol = &buffer[j*localHeight];
                T* col = &localBuffer[j*ldim];
                MemCopy( col, bufferCol, localHeight );
            }
        }
    }
};

template<typename T>
struct MakeUniformHelper<T,MD,STAR>
{
    static void Func( DistMatrix<T,MD,STAR>& A, T center, BASE(T) radius )
    {
        const Int n = A.Width();
        const Int localHeight = A.LocalHeight();
        for( Int j=0; j<n; ++j )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                A.SetLocal( iLoc, j, SampleBall( center, radius ) );
    }
};

template<typename T>
struct MakeUniformHelper<T,MR,MC>
{
    static void Func( DistMatrix<T,MR,MC>& A, T center, BASE(T) radius )
    {
        const Int localHeight = A.LocalHeight(); 
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                A.SetLocal( iLoc, jLoc, SampleBall( center, radius ) );
    }
};

template<typename T>
struct MakeUniformHelper<T,MR,STAR>
{
    static void Func( DistMatrix<T,MR,STAR>& A, T center, BASE(T) radius )
    {
        const Grid& grid = A.Grid();
        const Int n = A.Width();
        const Int localHeight = A.LocalHeight();
        const Int bufSize = localHeight*n;
        std::vector<T> buffer( bufSize );

        // Create random matrix on process row 0, then broadcast
        if( grid.Row() == 0 )
        {
            for( Int j=0; j<n; ++j )
                for( Int i=0; i<localHeight; ++i )
                    buffer[i+j*localHeight] = SampleBall( center, radius );
        }
        mpi::Broadcast( buffer.data(), bufSize, 0, grid.ColComm() );

        // Unpack
        T* localBuffer = A.Buffer();
        const Int ldim = A.LDim();
        PARALLEL_FOR COLLAPSE(2)
        for( Int j=0; j<n; ++j )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                localBuffer[iLoc+j*ldim] = buffer[iLoc+j*localHeight];
    }
};

template<typename T>
struct MakeUniformHelper<T,STAR,MC>
{
    static void Func( DistMatrix<T,STAR,MC>& A, T center, BASE(T) radius )
    {
        const Grid& grid = A.Grid();
        const Int m = A.Height();
        const Int localWidth = A.LocalWidth();
        const Int bufSize = m*localWidth;
        std::vector<T> buffer( bufSize );

        // Create a random matrix on process column 0, then broadcast
        if( grid.Col() == 0 )
        {
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                for( Int i=0; i<m; ++i )
                    buffer[i+jLoc*m] = SampleBall( center, radius );
        }
        mpi::Broadcast( buffer.data(), bufSize, 0, grid.RowComm() );

        // Unpack
        T* localBuffer = A.Buffer();
        const Int ldim = A.LDim();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const T* bufferCol = &buffer[jLoc*m];
            T* col = &localBuffer[jLoc*ldim];
            MemCopy( col, bufferCol, m );
        }
    }
};

template<typename T>
struct MakeUniformHelper<T,STAR,MD>
{
    static void Func( DistMatrix<T,STAR,MD>& A, T center, BASE(T) radius )
    {
        const Int m = A.Height();
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int i=0; i<m; ++i )
                A.SetLocal( i, jLoc, SampleBall( center, radius ) );
    }
};

template<typename T>
struct MakeUniformHelper<T,STAR,MR>
{
    static void Func( DistMatrix<T,STAR,MR>& A, T center, BASE(T) radius )
    {
        const Grid& grid = A.Grid();
        const Int m = A.Height();
        const Int localWidth = A.LocalWidth();
        const Int bufSize = m*localWidth;
        std::vector<T> buffer( bufSize );

        // Create random matrix on process row 0, then broadcast
        if( grid.Row() == 0 )
        {
            for( Int j=0; j<localWidth; ++j )
                for( Int i=0; i<m; ++i )
                    buffer[i+j*m] = SampleBall( center, radius );
        }
        mpi::Broadcast( buffer.data(), bufSize, 0, grid.ColComm() );

        // Unpack
        T* localBuffer = A.Buffer();
        const Int ldim = A.LDim();
        PARALLEL_FOR
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const T* bufferCol = &buffer[jLoc*m];
            T* col = &localBuffer[jLoc*ldim];
            MemCopy( col, bufferCol, m );
        }
    }
};

template<typename T>
struct MakeUniformHelper<T,STAR,STAR>
{
    static void Func( DistMatrix<T,STAR,STAR>& A, T center, BASE(T) radius )
    {
        const Grid& grid = A.Grid();
        const Int m = A.Height();
        const Int n = A.Width();
        const Int bufSize = m*n;

        if( grid.InGrid() )
        {
            std::vector<T> buffer( bufSize );

            if( grid.Rank() == 0 )
            {
                for( Int j=0; j<n; ++j )
                    for( Int i=0; i<m; ++i )
                        buffer[i+j*m] = SampleBall( center, radius );
            }
            mpi::Broadcast( buffer.data(), bufSize, 0, grid.Comm() );

            // Unpack
            T* localBuffer = A.Buffer();
            const Int ldim = A.LDim();
            PARALLEL_FOR
            for( Int j=0; j<n; ++j )
            {
                const T* bufferCol = &buffer[j*m];
                T* col = &localBuffer[j*ldim];
                MemCopy( col, bufferCol, m );
            }
        }
    }
};

template<typename T>
struct MakeUniformHelper<T,STAR,VC>
{
    static void Func( DistMatrix<T,STAR,VC>& A, T center, BASE(T) radius )
    {
        const Int m = A.Height();
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int i=0; i<m; ++i )
                A.SetLocal( i, jLoc, SampleBall( center, radius ) );
    }
};

template<typename T>
struct MakeUniformHelper<T,STAR,VR>
{
    static void Func( DistMatrix<T,STAR,VR>& A, T center, BASE(T) radius )
    {
        const Int m = A.Height();
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int i=0; i<m; ++i )
                A.SetLocal( i, jLoc, SampleBall( center, radius ) );
    }
};

template<typename T>
struct MakeUniformHelper<T,VC,STAR>
{
    static void Func( DistMatrix<T,VC,STAR>& A, T center, BASE(T) radius )
    {
        const Int n = A.Width();
        const Int localHeight = A.LocalHeight();
        for( Int j=0; j<n; ++j )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                A.SetLocal( iLoc, j, SampleBall( center, radius ) );
    }
};

template<typename T>
struct MakeUniformHelper<T,VR,STAR>
{
    static void Func( DistMatrix<T,VR,STAR>& A, T center, BASE(T) radius )
    {
        const Int n = A.Width();
        const Int localHeight = A.LocalHeight();
        for( Int j=0; j<n; ++j )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                A.SetLocal( iLoc, j, SampleBall( center, radius ) );
    }
};

} // namespace internal

template<typename T,Distribution U,Distribution V>
inline void
MakeUniform( DistMatrix<T,U,V>& A, T center=0, BASE(T) radius=1 )
{
#ifndef RELEASE
    CallStackEntry cse("Uniform");
#endif
    internal::MakeUniformHelper<T,U,V>::Func( A, center, radius );
}

template<typename T,Distribution U,Distribution V>
inline void
Uniform( DistMatrix<T,U,V>& A, Int m, Int n, T center=0, BASE(T) radius=1 )
{
#ifndef RELEASE
    CallStackEntry cse("Uniform");
#endif
    A.ResizeTo( m, n );
    MakeUniform( A, center, radius );
}

template<typename T,Distribution U=MC,Distribution V=MR>
inline DistMatrix<T,U,V>
Uniform( const Grid& g, Int m, Int n, T center=0, BASE(T) radius=1 )
{
    DistMatrix<T,U,V> A( m, n, g );
    MakeUniform( A, center, radius );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_UNIFORM_HPP
