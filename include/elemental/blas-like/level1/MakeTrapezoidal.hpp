/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_MAKETRAPEZOIDAL_HPP
#define ELEM_BLAS_MAKETRAPEZOIDAL_HPP

namespace elem {

template<typename T>
inline void
MakeTrapezoidal
( UpperOrLower uplo, Matrix<T>& A, Int offset=0, LeftOrRight side=LEFT )
{
#ifndef RELEASE
    CallStackEntry entry("MakeTrapezoidal");
#endif
    const Int height = A.Height();
    const Int width = A.Width();
    const Int ldim = A.LDim();
    T* buffer = A.Buffer();

    if( uplo == LOWER )
    {
        if( side == LEFT )
        {
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int j=Max(0,offset+1); j<width; ++j )
            {
                const Int lastZeroRow = j-offset-1;
                const Int numZeroRows = Min( lastZeroRow+1, height );
                MemZero( &buffer[j*ldim], numZeroRows );
            }
        }
        else
        {
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int j=Max(0,offset-height+width+1); j<width; ++j )
            {
                const Int lastZeroRow = j-offset+height-width-1;
                const Int numZeroRows = Min( lastZeroRow+1, height );
                MemZero( &buffer[j*ldim], numZeroRows );
            }
        }
    }
    else
    {
        if( side == LEFT )
        {
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int j=0; j<width; ++j )
            {
                const Int firstZeroRow = Max(j-offset+1,0);
                if( firstZeroRow < height )
                    MemZero
                    ( &buffer[firstZeroRow+j*ldim], height-firstZeroRow );
            }
        }
        else
        {
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for( Int j=0; j<width; ++j )
            {
                const Int firstZeroRow = Max(j-offset+height-width+1,0);
                if( firstZeroRow < height )
                    MemZero
                    ( &buffer[firstZeroRow+j*ldim], height-firstZeroRow );
            }
        }
    }
}

template<typename T,Distribution U,Distribution V>
inline void
MakeTrapezoidal
( UpperOrLower uplo, DistMatrix<T,U,V>& A, Int offset=0, LeftOrRight side=LEFT )
{
#ifndef RELEASE
    CallStackEntry entry("MakeTrapezoidal");
#endif
    const Int height = A.Height();
    const Int width = A.Width();
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int colShift = A.ColShift();
    const Int rowShift = A.RowShift();
    const Int colStride = A.ColStride();
    const Int rowStride = A.RowStride();

    T* buffer = A.Buffer();
    const Int ldim = A.LDim();

    if( uplo == LOWER )
    {

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = rowShift + jLoc*rowStride;
            const Int lastZeroRow =
                ( side==LEFT ? j-offset-1
                             : j-offset+height-width-1 );
            if( lastZeroRow >= 0 )
            {
                const Int boundary = Min( lastZeroRow+1, height );
                const Int numZeroRows =
                    Length_( boundary, colShift, colStride );
                MemZero( &buffer[jLoc*ldim], numZeroRows );
            }
        }
    }
    else
    {
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = rowShift + jLoc*rowStride;
            const Int firstZeroRow =
                ( side==LEFT ? Max(j-offset+1,0)
                             : Max(j-offset+height-width+1,0) );
            const Int numNonzeroRows =
                Length_(firstZeroRow,colShift,colStride);
            if( numNonzeroRows < localHeight )
            {
                T* col = &buffer[numNonzeroRows+jLoc*ldim];
                MemZero( col, localHeight-numNonzeroRows );
            }
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_BLAS_MAKETRAPEZOIDAL_HPP
