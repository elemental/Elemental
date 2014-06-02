/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_GAUSSIAN_HPP
#define EL_GAUSSIAN_HPP

namespace El {

// Draw each entry from a normal PDF
template<typename T>
inline void
MakeGaussian( Matrix<T>& A, T mean=0, Base<T> stddev=1 )
{
    DEBUG_ONLY(CallStackEntry cse("MakeGaussian"))
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, SampleNormal( mean, stddev ) );
}

template<typename T>
inline void
MakeGaussian( AbstractDistMatrix<T>& A, T mean=0, Base<T> stddev=1 )
{
    DEBUG_ONLY(CallStackEntry cse("MakeGaussian"))
    if( A.RedundantSize() == 1 )
    {
        MakeGaussian( A.Matrix(), mean, stddev );
    }
    else if( A.Participating() && A.LocalHeight() == A.LDim() )
    {
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        if( A.RedundantRank() == 0 )
            MakeGaussian( A.Matrix(), mean, stddev );
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
                        SampleNormal( mean, stddev );
        }
        mpi::Broadcast( buffer.data(), bufSize, 0, A.RedundantComm() );
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( A.Buffer(0,jLoc), &buffer[jLoc*localHeight], localHeight );
    }
}

template<typename T>
inline void
MakeGaussian( AbstractBlockDistMatrix<T>& A, T mean=0, Base<T> stddev=1 )
{
    DEBUG_ONLY(CallStackEntry cse("MakeGaussian"))
    if( A.RedundantSize() == 1 )
    {
        MakeGaussian( A.Matrix(), mean, stddev );
    }
    else if( A.Participating() && A.LocalHeight() == A.LDim() )
    {
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        if( A.RedundantRank() == 0 )
            MakeGaussian( A.Matrix(), mean, stddev );
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
                        SampleNormal( mean, stddev );
        }
        mpi::Broadcast( buffer.data(), bufSize, 0, A.RedundantComm() );
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy
            ( A.Buffer(0,jLoc), &buffer[jLoc*localHeight], localHeight );
    }
}

template<typename T>
inline void
Gaussian( Matrix<T>& A, Int m, Int n, T mean=0, Base<T> stddev=1 )
{
    DEBUG_ONLY(CallStackEntry cse("Gaussian"))
    A.Resize( m, n );
    MakeGaussian( A, mean, stddev );
}

template<typename T>
inline void
Gaussian( AbstractDistMatrix<T>& A, Int m, Int n, T mean=0, Base<T> stddev=1 )
{
    DEBUG_ONLY(CallStackEntry cse("Gaussian"))
    A.Resize( m, n );
    MakeGaussian( A, mean, stddev );
}

template<typename T>
inline void
Gaussian
( AbstractBlockDistMatrix<T>& A, Int m, Int n, T mean=0, Base<T> stddev=1 )
{
    DEBUG_ONLY(CallStackEntry cse("Gaussian"))
    A.Resize( m, n );
    MakeGaussian( A, mean, stddev );
}

} // namespace El

#endif // ifndef EL_GAUSSIAN_HPP
