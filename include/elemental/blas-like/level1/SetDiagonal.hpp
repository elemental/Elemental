/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_SETDIAGONAL_HPP
#define ELEM_BLAS_SETDIAGONAL_HPP

namespace elem {

template<typename T>
inline void
SetDiagonal( Matrix<T>& A, T alpha )
{
#ifndef RELEASE
    CallStackEntry entry("SetDiagonal");
#endif
    const int height = A.Height();
    const int width = A.Width();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( int j=0; j<std::min(height,width); ++j )
        A.Set(j,j,alpha);
}

template<typename T>
inline void
SetDiagonal( Matrix<T>& A, T alpha, int offset, LeftOrRight side=LEFT )
{
#ifndef RELEASE
    CallStackEntry entry("SetDiagonal");
#endif
    const int height = A.Height();
    const int width = A.Width();
    if( side == LEFT )
    {
        for( int j=0; j<width; ++j )
        {
            const int i = j-offset;
            if( i >= 0 && i < height )
                A.Set(i,j,alpha);
        }
    }
    else
    {
        for( int j=0; j<width; ++j )
        {
            const int i = j-offset+height-width;
            if( i >= 0 && i < height )
                A.Set(i,j,alpha);
        }
    }
}

template<typename T,Distribution U,Distribution V>
inline void
SetDiagonal( DistMatrix<T,U,V>& A, T alpha )
{
#ifndef RELEASE
    CallStackEntry entry("SetDiagonal");
#endif
    const int height = A.Height();
    const int rowShift = A.RowShift();
    const int colShift = A.ColShift();
    const int rowStride = A.RowStride();
    const int colStride = A.ColStride();
    const int localWidth = A.LocalWidth();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        if( j < height && j % colStride == colShift )
        {
            const int iLoc = (j-colShift) / colStride;
            A.SetLocal( iLoc, jLoc, alpha );
        }
    }
}

template<typename T,Distribution U,Distribution V>
inline void
SetDiagonal
( DistMatrix<T,U,V>& A, T alpha, int offset, LeftOrRight side=LEFT )
{
#ifndef RELEASE
    CallStackEntry entry("SetDiagonal");
#endif
    const int height = A.Height();
    const int width = A.Width();
    const int rowShift = A.RowShift();
    const int colShift = A.ColShift();
    const int rowStride = A.RowStride();
    const int colStride = A.ColStride();
    const int localWidth = A.LocalWidth();
    if( side == LEFT )
    {
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*rowStride;
            const int i = j-offset;
            if( i >= 0 && i < height && i % colStride == colShift )
            {
                const int iLoc = (i-colShift) / colStride;
                A.SetLocal( iLoc, jLoc, alpha );
            }
        }
    }
    else
    {
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*rowStride;
            const int i = j-offset+height-width;
            if( i >= 0 && i < height && i % colStride == colShift )
            {
                const int iLoc = (i-colShift) / colStride;
                A.SetLocal( iLoc, jLoc, alpha );
            }
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_BLAS_SETDIAGONAL_HPP
