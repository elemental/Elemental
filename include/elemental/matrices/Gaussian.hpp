/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_GAUSSIAN_HPP
#define ELEM_MATRICES_GAUSSIAN_HPP

namespace elem {

// Draw each entry from a normal PDF
template<typename T>
inline void
MakeGaussian( Matrix<T>& A, T mean=0, BASE(T) stddev=1 )
{
#ifndef RELEASE
    CallStackEntry cse("MakeGaussian");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, Normal( mean, stddev ) );
}

template<typename T>
inline void
Gaussian( Matrix<T>& A, Int m, Int n, T mean=0, BASE(T) stddev=1 )
{
#ifndef RELEASE
    CallStackEntry cse("Gaussian");
#endif
    A.ResizeTo( m, n );
    MakeGaussian( A, mean, stddev );
}

template<typename T>
inline Matrix<T>
Gaussian( Int m, Int n, T mean=0, BASE(T) stddev=1 )
{
    Matrix<T> A( m, n );
    MakeGaussian( A, mean, stddev );
    return A;
}

namespace internal {

template<typename T,Distribution U,Distribution V>
struct MakeGaussianHelper
{
    static void Func( DistMatrix<T,U,V>& A, T mean, BASE(T) stddev );  
};

template<typename T>
struct MakeGaussianHelper<T,CIRC,CIRC>
{
    static void Func( DistMatrix<T,CIRC,CIRC>& A, T mean, BASE(T) stddev )
    {
        if( A.Grid().VCRank() == A.Root() )
        {
            const Int height = A.Height(); 
            const Int width = A.Width();
            for( Int j=0; j<width; ++j )
                for( Int i=0; i<height; ++i )
                    A.SetLocal( i, j, Normal( mean, stddev ) );
        }
    }
};

template<typename T>
struct MakeGaussianHelper<T,MC,MR>
{
    static void Func( DistMatrix<T,MC,MR>& A, T mean, BASE(T) stddev )
    {
        const Int localHeight = A.LocalHeight(); 
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                A.SetLocal( iLoc, jLoc, Normal( mean, stddev ) );
    }
};

template<typename T>
struct MakeGaussianHelper<T,MC,STAR>
{
    static void Func( DistMatrix<T,MC,STAR>& A, T mean, BASE(T) stddev )
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
                        buffer[iLoc+j*localHeight] = Normal( mean, stddev );
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
struct MakeGaussianHelper<T,MD,STAR>
{
    static void Func( DistMatrix<T,MD,STAR>& A, T mean, BASE(T) stddev )
    {
        const Int n = A.Width();
        const Int localHeight = A.LocalHeight();
        for( Int j=0; j<n; ++j )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                A.SetLocal( iLoc, j, Normal( mean, stddev ) );
    }
};

template<typename T>
struct MakeGaussianHelper<T,MR,MC>
{
    static void Func( DistMatrix<T,MR,MC>& A, T mean, BASE(T) stddev )
    {
        const Int localHeight = A.LocalHeight(); 
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                A.SetLocal( iLoc, jLoc, Normal( mean, stddev ) );
    }
};

template<typename T>
struct MakeGaussianHelper<T,MR,STAR>
{
    static void Func( DistMatrix<T,MR,STAR>& A, T mean, BASE(T) stddev )
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
                    buffer[i+j*localHeight] = Normal( mean, stddev );
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
struct MakeGaussianHelper<T,STAR,MC>
{
    static void Func( DistMatrix<T,STAR,MC>& A, T mean, BASE(T) stddev )
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
                    buffer[i+jLoc*m] = Normal( mean, stddev );
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
struct MakeGaussianHelper<T,STAR,MD>
{
    static void Func( DistMatrix<T,STAR,MD>& A, T mean, BASE(T) stddev )
    {
        const Int m = A.Height();
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int i=0; i<m; ++i )
                A.SetLocal( i, jLoc, Normal( mean, stddev ) );
    }
};

template<typename T>
struct MakeGaussianHelper<T,STAR,MR>
{
    static void Func( DistMatrix<T,STAR,MR>& A, T mean, BASE(T) stddev )
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
                    buffer[i+j*m] = Normal( mean, stddev );
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
struct MakeGaussianHelper<T,STAR,STAR>
{
    static void Func( DistMatrix<T,STAR,STAR>& A, T mean, BASE(T) stddev )
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
                        buffer[i+j*m] = Normal( mean, stddev );
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
struct MakeGaussianHelper<T,STAR,VC>
{
    static void Func( DistMatrix<T,STAR,VC>& A, T mean, BASE(T) stddev )
    {
        const Int m = A.Height();
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int i=0; i<m; ++i )
                A.SetLocal( i, jLoc, Normal( mean, stddev ) );
    }
};

template<typename T>
struct MakeGaussianHelper<T,STAR,VR>
{
    static void Func( DistMatrix<T,STAR,VR>& A, T mean, BASE(T) stddev )
    {
        const Int m = A.Height();
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            for( Int i=0; i<m; ++i )
                A.SetLocal( i, jLoc, Normal( mean, stddev ) );
    }
};

template<typename T>
struct MakeGaussianHelper<T,VC,STAR>
{
    static void Func( DistMatrix<T,VC,STAR>& A, T mean, BASE(T) stddev )
    {
        const Int n = A.Width();
        const Int localHeight = A.LocalHeight();
        for( Int j=0; j<n; ++j )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                A.SetLocal( iLoc, j, Normal( mean, stddev ) );
    }
};

template<typename T>
struct MakeGaussianHelper<T,VR,STAR>
{
    static void Func( DistMatrix<T,VR,STAR>& A, T mean, BASE(T) stddev )
    {
        const Int n = A.Width();
        const Int localHeight = A.LocalHeight();
        for( Int j=0; j<n; ++j )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                A.SetLocal( iLoc, j, Normal( mean, stddev ) );
    }
};

} // namespace internal

template<typename T,Distribution U,Distribution V>
inline void
MakeGaussian( DistMatrix<T,U,V>& A, T mean=0, BASE(T) stddev=1 )
{
#ifndef RELEASE
    CallStackEntry cse("Gaussian");
#endif
    internal::MakeGaussianHelper<T,U,V>::Func( A, mean, stddev );
}

template<typename T,Distribution U,Distribution V>
inline void
Gaussian( DistMatrix<T,U,V>& A, Int m, Int n, T mean=0, BASE(T) stddev=1 )
{
#ifndef RELEASE
    CallStackEntry cse("Gaussian");
#endif
    A.ResizeTo( m, n );
    MakeGaussian( A, mean, stddev );
}

template<typename T,Distribution U=MC,Distribution V=MR>
inline DistMatrix<T,U,V>
Gaussian( const Grid& g, Int m, Int n, T mean=0, BASE(T) stddev=1 )
{
    DistMatrix<T,U,V> A( m, n, g );
    MakeGaussian( A, mean, stddev );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_GAUSSIAN_HPP
