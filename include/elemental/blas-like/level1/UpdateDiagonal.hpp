/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_UPDATEDIAGONAL_HPP
#define ELEM_BLAS_UPDATEDIAGONAL_HPP

// This is essentially equivalent to SetDiagonal, but with s/Set/Update/g.

namespace elem {

template<typename T>
inline void
UpdateDiagonal( Matrix<T>& A, T alpha )
{
#ifndef RELEASE
    CallStackEntry entry("UpdateDiagonal");
#endif
    const Int height = A.Height();
    const Int width = A.Width();
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for( Int j=0; j<Min(height,width); ++j )
        A.Update(j,j,alpha);
}

template<typename T>
inline void
UpdateDiagonal( Matrix<T>& A, T alpha, Int offset, LeftOrRight side=LEFT )
{
#ifndef RELEASE
    CallStackEntry entry("UpdateDiagonal");
#endif
    const Int height = A.Height();
    const Int width = A.Width();
    if( side == LEFT )
    {
        for( Int j=0; j<width; ++j )
        {
            const Int i = j-offset;
            if( i >= 0 && i < height )
                A.Update(i,j,alpha);
        }
    }
    else
    {
        for( Int j=0; j<width; ++j )
        {
            const Int i = j-offset+height-width;
            if( i >= 0 && i < height )
                A.Update(i,j,alpha);
        }
    }
}

template<typename T,Distribution U,Distribution V>
inline void
UpdateDiagonal( DistMatrix<T,U,V>& A, T alpha )
{
#ifndef RELEASE
    CallStackEntry entry("UpdateDiagonal");
#endif
    const Int height = A.Height();
    const Int rowShift = A.RowShift();
    const Int colShift = A.ColShift();
    const Int rowStride = A.RowStride();
    const Int colStride = A.ColStride();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        if( j < height && j % colStride == colShift )
        {
            const Int iLoc = (j-colShift) / colStride;
            A.UpdateLocal( iLoc, jLoc, alpha );
        }
    }
}

template<typename T,Distribution U,Distribution V>
inline void
UpdateDiagonal
( DistMatrix<T,U,V>& A, T alpha, Int offset, LeftOrRight side=LEFT )
{
#ifndef RELEASE
    CallStackEntry entry("UpdateDiagonal");
#endif
    const Int height = A.Height();
    const Int width = A.Width();
    const Int rowShift = A.RowShift();
    const Int colShift = A.ColShift();
    const Int rowStride = A.RowStride();
    const Int colStride = A.ColStride();
    const Int localWidth = A.LocalWidth();
    if( side == LEFT )
    {
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = rowShift + jLoc*rowStride;
            const Int i = j-offset;
            if( i >= 0 && i < height && i % colStride == colShift )
            {
                const Int iLoc = (i-colShift) / colStride;
                A.UpdateLocal( iLoc, jLoc, alpha );
            }
        }
    }
    else
    {
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = rowShift + jLoc*rowStride;
            const Int i = j-offset+height-width;
            if( i >= 0 && i < height && i % colStride == colShift )
            {
                const Int iLoc = (i-colShift) / colStride;
                A.UpdateLocal( iLoc, jLoc, alpha );
            }
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_BLAS_UPDATEDIAGONAL_HPP
