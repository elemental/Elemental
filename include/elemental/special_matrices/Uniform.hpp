/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

// Draw each entry from a uniform PDF over the closed unit ball.
template<typename T>
inline void
MakeUniform( Matrix<T>& A, T center, typename Base<T>::type radius )
{
#ifndef RELEASE
    PushCallStack("MakeUniform");
#endif
    const int m = A.Height();
    const int n = A.Width();
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            A.Set( i, j, center+radius*SampleUnitBall<T>() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Uniform( int m, int n, Matrix<T>& A, T center, typename Base<T>::type radius )
{
#ifndef RELEASE
    PushCallStack("Uniform");
#endif
    A.ResizeTo( m, n );
    MakeUniform( A, center, radius );
#ifndef RELEASE
    PopCallStack();
#endif
}

namespace internal {

template<typename T,Distribution U,Distribution V>
struct MakeUniformHelper
{
    static void Func
    ( DistMatrix<T,U,V>& A, T center, typename Base<T>::type radius );  
};

template<typename T>
struct MakeUniformHelper<T,MC,MR>
{
    static void Func
    ( DistMatrix<T,MC,MR>& A, T center, typename Base<T>::type radius )
    {
        const int localHeight = A.LocalHeight(); 
        const int localWidth = A.LocalWidth();
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                A.SetLocal( iLocal, jLocal, center+radius*SampleUnitBall<T>() );
    }
};

template<typename T>
struct MakeUniformHelper<T,MC,STAR>
{
    static void Func
    ( DistMatrix<T,MC,STAR>& A, T center, typename Base<T>::type radius )
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
                    for( int iLocal=0; iLocal<localHeight; ++iLocal )
                        buffer[iLocal+j*localHeight] = 
                            center + radius*SampleUnitBall<T>();
            }
            mpi::Broadcast( &buffer[0], bufSize, 0, grid.RowComm() );

            // Unpack
            T* localBuffer = A.LocalBuffer();
            const int ldim = A.LocalLDim();
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
    static void Func
    ( DistMatrix<T,MD,STAR>& A, T center, typename Base<T>::type radius )
    {
        if( A.Participating() )
        {
            const int n = A.Width();
            const int localHeight = A.LocalHeight();
            for( int j=0; j<n; ++j )
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                    A.SetLocal( iLocal, j, center+radius*SampleUnitBall<T>() );
        }
    }
};

template<typename T>
struct MakeUniformHelper<T,MR,MC>
{
    static void Func
    ( DistMatrix<T,MR,MC>& A, T center, typename Base<T>::type radius )
    {
        const int localHeight = A.LocalHeight(); 
        const int localWidth = A.LocalWidth();
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                A.SetLocal( iLocal, jLocal, center+radius*SampleUnitBall<T>() );
    }
};

template<typename T>
struct MakeUniformHelper<T,MR,STAR>
{
    static void Func
    ( DistMatrix<T,MR,STAR>& A, T center, typename Base<T>::type radius )
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
        T* localBuffer = A.LocalBuffer();
        const int ldim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for COLLAPSE(2)
#endif
        for( int j=0; j<n; ++j )
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                localBuffer[iLocal+j*ldim] = buffer[iLocal+j*localHeight];
    }
};

template<typename T>
struct MakeUniformHelper<T,STAR,MC>
{
    static void Func
    ( DistMatrix<T,STAR,MC>& A, T center, typename Base<T>::type radius )
    {
        const Grid& grid = A.Grid();
        const int m = A.Height();
        const int localWidth = A.LocalWidth();
        const int bufSize = m*localWidth;
        std::vector<T> buffer( bufSize );

        // Create a random matrix on process column 0, then broadcast
        if( grid.Col() == 0 )
        {
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                for( int i=0; i<m; ++i )
                    buffer[i+jLocal*m] = center+radius*SampleUnitBall<T>();
        }
        mpi::Broadcast( &buffer[0], bufSize, 0, grid.RowComm() );

        // Unpack
        T* localBuffer = A.LocalBuffer();
        const int ldim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* bufferCol = &buffer[jLocal*m];
            T* col = &localBuffer[jLocal*ldim];
            MemCopy( col, bufferCol, m );
        }
    }
};

template<typename T>
struct MakeUniformHelper<T,STAR,MD>
{
    static void Func
    ( DistMatrix<T,STAR,MD>& A, T center, typename Base<T>::type radius )
    {
        if( A.Participating() )
        {
            const int m = A.Height();
            const int localWidth = A.LocalWidth();
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                for( int i=0; i<m; ++i )
                    A.SetLocal( i, jLocal, center+radius*SampleUnitBall<T>() );
        }
    }
};

template<typename T>
struct MakeUniformHelper<T,STAR,MR>
{
    static void Func
    ( DistMatrix<T,STAR,MR>& A, T center, typename Base<T>::type radius )
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
        T* localBuffer = A.LocalBuffer();
        const int ldim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const T* bufferCol = &buffer[jLocal*m];
            T* col = &localBuffer[jLocal*ldim];
            MemCopy( col, bufferCol, m );
        }
    }
};

template<typename T>
struct MakeUniformHelper<T,STAR,STAR>
{
    static void Func
    ( DistMatrix<T,STAR,STAR>& A, T center, typename Base<T>::type radius )
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
            T* localBuffer = A.LocalBuffer();
            const int ldim = A.LocalLDim();
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
    static void Func
    ( DistMatrix<T,STAR,VC>& A, T center, typename Base<T>::type radius )
    {
        const int m = A.Height();
        const int localWidth = A.LocalWidth();
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
            for( int i=0; i<m; ++i )
                A.SetLocal( i, jLocal, center+radius*SampleUnitBall<T>() );
    }
};

template<typename T>
struct MakeUniformHelper<T,STAR,VR>
{
    static void Func
    ( DistMatrix<T,STAR,VR>& A, T center, typename Base<T>::type radius )
    {
        const int m = A.Height();
        const int localWidth = A.LocalWidth();
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
            for( int i=0; i<m; ++i )
                A.SetLocal( i, jLocal, center+radius*SampleUnitBall<T>() );
    }
};

template<typename T>
struct MakeUniformHelper<T,VC,STAR>
{
    static void Func
    ( DistMatrix<T,VC,STAR>& A, T center, typename Base<T>::type radius )
    {
        const int n = A.Width();
        const int localHeight = A.LocalHeight();
        for( int j=0; j<n; ++j )
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                A.SetLocal( iLocal, j, center+radius*SampleUnitBall<T>() );
    }
};

template<typename T>
struct MakeUniformHelper<T,VR,STAR>
{
    static void Func
    ( DistMatrix<T,VR,STAR>& A, T center, typename Base<T>::type radius )
    {
        const int n = A.Width();
        const int localHeight = A.LocalHeight();
        for( int j=0; j<n; ++j )
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
                A.SetLocal( iLocal, j, center+radius*SampleUnitBall<T>() );
    }
};

} // namespace internal

template<typename T,Distribution U,Distribution V>
inline void
MakeUniform
( DistMatrix<T,U,V>& A, T center, typename Base<T>::type radius )
{
#ifndef RELEASE
    PushCallStack("Uniform");
#endif
    internal::MakeUniformHelper<T,U,V>::Func( A, center, radius );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
Uniform
( int m, int n, DistMatrix<T,U,V>& A, T center, typename Base<T>::type radius )
{
#ifndef RELEASE
    PushCallStack("Uniform");
#endif
    A.ResizeTo( m, n );
    MakeUniform( A, center, radius );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
