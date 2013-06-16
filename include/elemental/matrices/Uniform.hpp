/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_UNIFORM_HPP
#define MATRICES_UNIFORM_HPP

namespace elem {

// Draw each entry from a uniform PDF over the closed unit ball.
template<typename T>
inline void
MakeUniform( Matrix<T>& A, T center=0, BASE(T) radius=1 )
{
#ifndef RELEASE
    CallStackEntry entry("MakeUniform");
#endif
    const int m = A.Height();
    const int n = A.Width();
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            A.Set( i, j, center+radius*SampleUnitBall<T>() );
}

template<typename T>
inline void
Uniform
( Matrix<T>& A, int m, int n, T center=0, BASE(T) radius=1 )
{
#ifndef RELEASE
    CallStackEntry entry("Uniform");
#endif
    A.ResizeTo( m, n );
    MakeUniform( A, center, radius );
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
            const int height = A.Height(); 
            const int width = A.Width();
            for( int j=0; j<width; ++j )
                for( int i=0; i<height; ++i )
                    A.SetLocal( i, j, center+radius*SampleUnitBall<T>() );
        }
    }
};

template<typename T>
struct MakeUniformHelper<T,MC,MR>
{
    static void Func( DistMatrix<T,MC,MR>& A, T center, BASE(T) radius )
    {
        const int localHeight = A.LocalHeight(); 
        const int localWidth = A.LocalWidth();
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
            for( int iLoc=0; iLoc<localHeight; ++iLoc )
                A.SetLocal( iLoc, jLoc, center+radius*SampleUnitBall<T>() );
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
            const int n = A.Width();
            const int localHeight = A.LocalHeight();
            const int bufSize = localHeight*n;
            std::vector<T> buffer( bufSize );

            // Create random matrix on process column 0, then broadcast
            if( grid.Col() == 0 )
            {
                for( int j=0; j<n; ++j )
                    for( int iLoc=0; iLoc<localHeight; ++iLoc )
                        buffer[iLoc+j*localHeight] = 
                            center + radius*SampleUnitBall<T>();
            }
            mpi::Broadcast( &buffer[0], bufSize, 0, grid.RowComm() );

            // Unpack
            T* localBuffer = A.Buffer();
            const int ldim = A.LDim();
#ifdef HAVE_OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<n; ++j )
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
        if( A.Participating() )
        {
            const int n = A.Width();
            const int localHeight = A.LocalHeight();
            for( int j=0; j<n; ++j )
                for( int iLoc=0; iLoc<localHeight; ++iLoc )
                    A.SetLocal( iLoc, j, center+radius*SampleUnitBall<T>() );
        }
    }
};

template<typename T>
struct MakeUniformHelper<T,MR,MC>
{
    static void Func( DistMatrix<T,MR,MC>& A, T center, BASE(T) radius )
    {
        const int localHeight = A.LocalHeight(); 
        const int localWidth = A.LocalWidth();
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
            for( int iLoc=0; iLoc<localHeight; ++iLoc )
                A.SetLocal( iLoc, jLoc, center+radius*SampleUnitBall<T>() );
    }
};

template<typename T>
struct MakeUniformHelper<T,MR,STAR>
{
    static void Func( DistMatrix<T,MR,STAR>& A, T center, BASE(T) radius )
    {
        const Grid& grid = A.Grid();
        const int n = A.Width();
        const int localHeight = A.LocalHeight();
        const int bufSize = localHeight*n;
        std::vector<T> buffer( bufSize );

        // Create random matrix on process row 0, then broadcast
        if( grid.Row() == 0 )
        {
            for( int j=0; j<n; ++j )
                for( int i=0; i<localHeight; ++i )
                    buffer[i+j*localHeight] = center+radius*SampleUnitBall<T>();
        }
        mpi::Broadcast( &buffer[0], bufSize, 0, grid.ColComm() );

        // Unpack
        T* localBuffer = A.Buffer();
        const int ldim = A.LDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<n; ++j )
            for( int iLoc=0; iLoc<localHeight; ++iLoc )
                localBuffer[iLoc+j*ldim] = buffer[iLoc+j*localHeight];
    }
};

template<typename T>
struct MakeUniformHelper<T,STAR,MC>
{
    static void Func( DistMatrix<T,STAR,MC>& A, T center, BASE(T) radius )
    {
        const Grid& grid = A.Grid();
        const int m = A.Height();
        const int localWidth = A.LocalWidth();
        const int bufSize = m*localWidth;
        std::vector<T> buffer( bufSize );

        // Create a random matrix on process column 0, then broadcast
        if( grid.Col() == 0 )
        {
            for( int jLoc=0; jLoc<localWidth; ++jLoc )
                for( int i=0; i<m; ++i )
                    buffer[i+jLoc*m] = center+radius*SampleUnitBall<T>();
        }
        mpi::Broadcast( &buffer[0], bufSize, 0, grid.RowComm() );

        // Unpack
        T* localBuffer = A.Buffer();
        const int ldim = A.LDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
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
        if( A.Participating() )
        {
            const int m = A.Height();
            const int localWidth = A.LocalWidth();
            for( int jLoc=0; jLoc<localWidth; ++jLoc )
                for( int i=0; i<m; ++i )
                    A.SetLocal( i, jLoc, center+radius*SampleUnitBall<T>() );
        }
    }
};

template<typename T>
struct MakeUniformHelper<T,STAR,MR>
{
    static void Func( DistMatrix<T,STAR,MR>& A, T center, BASE(T) radius )
    {
        const Grid& grid = A.Grid();
        const int m = A.Height();
        const int localWidth = A.LocalWidth();
        const int bufSize = m*localWidth;
        std::vector<T> buffer( bufSize );

        // Create random matrix on process row 0, then broadcast
        if( grid.Row() == 0 )
        {
            for( int j=0; j<localWidth; ++j )
                for( int i=0; i<m; ++i )
                    buffer[i+j*m] = center+radius*SampleUnitBall<T>();
        }
        mpi::Broadcast( &buffer[0], bufSize, 0, grid.ColComm() );

        // Unpack
        T* localBuffer = A.Buffer();
        const int ldim = A.LDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
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
        const int m = A.Height();
        const int n = A.Width();
        const int bufSize = m*n;

        if( grid.InGrid() )
        {
            std::vector<T> buffer( bufSize );

            if( grid.Rank() == 0 )
            {
                for( int j=0; j<n; ++j )
                    for( int i=0; i<m; ++i )
                        buffer[i+j*m] = center+radius*SampleUnitBall<T>();
            }
            mpi::Broadcast( &buffer[0], bufSize, 0, grid.Comm() );

            // Unpack
            T* localBuffer = A.Buffer();
            const int ldim = A.LDim();
#ifdef HAVE_OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<n; ++j )
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
        const int m = A.Height();
        const int localWidth = A.LocalWidth();
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
            for( int i=0; i<m; ++i )
                A.SetLocal( i, jLoc, center+radius*SampleUnitBall<T>() );
    }
};

template<typename T>
struct MakeUniformHelper<T,STAR,VR>
{
    static void Func( DistMatrix<T,STAR,VR>& A, T center, BASE(T) radius )
    {
        const int m = A.Height();
        const int localWidth = A.LocalWidth();
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
            for( int i=0; i<m; ++i )
                A.SetLocal( i, jLoc, center+radius*SampleUnitBall<T>() );
    }
};

template<typename T>
struct MakeUniformHelper<T,VC,STAR>
{
    static void Func( DistMatrix<T,VC,STAR>& A, T center, BASE(T) radius )
    {
        const int n = A.Width();
        const int localHeight = A.LocalHeight();
        for( int j=0; j<n; ++j )
            for( int iLoc=0; iLoc<localHeight; ++iLoc )
                A.SetLocal( iLoc, j, center+radius*SampleUnitBall<T>() );
    }
};

template<typename T>
struct MakeUniformHelper<T,VR,STAR>
{
    static void Func( DistMatrix<T,VR,STAR>& A, T center, BASE(T) radius )
    {
        const int n = A.Width();
        const int localHeight = A.LocalHeight();
        for( int j=0; j<n; ++j )
            for( int iLoc=0; iLoc<localHeight; ++iLoc )
                A.SetLocal( iLoc, j, center+radius*SampleUnitBall<T>() );
    }
};

} // namespace internal

template<typename T,Distribution U,Distribution V>
inline void
MakeUniform( DistMatrix<T,U,V>& A, T center=0, BASE(T) radius=1 )
{
#ifndef RELEASE
    CallStackEntry entry("Uniform");
#endif
    internal::MakeUniformHelper<T,U,V>::Func( A, center, radius );
}

template<typename T,Distribution U,Distribution V>
inline void
Uniform( DistMatrix<T,U,V>& A, int m, int n, T center=0, BASE(T) radius=1 )
{
#ifndef RELEASE
    CallStackEntry entry("Uniform");
#endif
    A.ResizeTo( m, n );
    MakeUniform( A, center, radius );
}

} // namespace elem

#endif // ifndef MATRICES_UNIFORM_HPP
