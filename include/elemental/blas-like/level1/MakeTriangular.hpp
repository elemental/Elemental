/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_MAKETRIANGULAR_HPP
#define BLAS_MAKETRIANGULAR_HPP

namespace elem {

template<typename T>
inline void
MakeTriangular( UpperOrLower uplo, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("MakeTriangular");
#endif
    const int height = A.Height();
    const int width = A.Width();
    const int ldim = A.LDim();
    T* buffer = A.Buffer();

    if( uplo == LOWER )
    {
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( int j=1; j<width; ++j )
        {
            const int numZeroRows = std::min( j, height );
            MemZero( &buffer[j*ldim], numZeroRows );
        }
    }
    else
    {
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( int j=0; j<std::min(width,height); ++j )
        {
            const int firstZeroRow = j+1;
            MemZero( &buffer[firstZeroRow+j*ldim], height-firstZeroRow );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
MakeTriangular( UpperOrLower uplo, DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("MakeTriangular");
#endif
    const int height = A.Height();
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();
    const int colStride = A.ColStride();
    const int rowStride = A.RowStride();

    T* localBuffer = A.LocalBuffer();
    const int ldim = A.LocalLDim();

    if( uplo == LOWER )
    {

#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const int j = rowShift + jLocal*rowStride;
            const int lastZeroRow = j-1;
            if( lastZeroRow >= 0 )
            {
                const int boundary = std::min( lastZeroRow+1, height );
                const int numZeroRows =
                    RawLocalLength( boundary, colShift, colStride );
                MemZero( &localBuffer[jLocal*ldim], numZeroRows );
            }
        }
    }
    else
    {
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            const int j = rowShift + jLocal*rowStride;
            const int firstZeroRow = j+1;
            const int numNonzeroRows =
                RawLocalLength(firstZeroRow,colShift,colStride);
            if( numNonzeroRows < localHeight )
            {
                T* col = &localBuffer[numNonzeroRows+jLocal*ldim];
                MemZero( col, localHeight-numNonzeroRows );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef BLAS_MAKETRIANGULAR_HPP
