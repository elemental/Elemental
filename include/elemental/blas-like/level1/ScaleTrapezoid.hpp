/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

template<typename T>
inline void
ScaleTrapezoid
( T alpha, LeftOrRight side, UpperOrLower uplo, int offset, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("ScaleTrapezoid");
#endif
    const int height = A.Height();
    const int width = A.Width();
    const int ldim = A.LDim();
    T* buffer = A.Buffer();

    if( uplo == UPPER )
    {
        if( side == LEFT )
        {
#ifdef HAVE_OPENMP
            #pragma omp parallel for
#endif
            for( int j=std::max(0,offset-1); j<width; ++j )
            {
                const int numRows = j-offset+1;
                for( int i=0; i<numRows; ++i )
                    buffer[i+j*ldim] *= alpha;
            }
        }
        else
        {
#ifdef HAVE_OPENMP
            #pragma omp parallel for
#endif
            for( int j=std::max(0,offset-height+width-1); j<width; ++j )
            {
                const int numRows = j-offset+height-width+1;
                for( int i=0; i<numRows; ++i )
                    buffer[i+j*ldim] *= alpha;
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
            for( int j=0; j<width; ++j )
            {
                const int numZeroRows = std::max(j-offset,0);
                for( int i=numZeroRows; i<height; ++i )
                    buffer[i+j*ldim] *= alpha;
            }
        }
        else
        {
#ifdef HAVE_OPENMP
            #pragma omp parallel for
#endif
            for( int j=0; j<width; ++j )
            {
                const int numZeroRows = std::max(j-offset+height-width,0);
                for( int i=numZeroRows; i<height; ++i )
                    buffer[i+j*ldim] *= alpha;
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
ScaleTrapezoid
( T alpha, LeftOrRight side, UpperOrLower uplo, int offset,
  DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("ScaleTrapezoid");
#endif
    const int height = A.Height();
    const int width = A.Width();
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();
    const int colStride = A.ColStride();
    const int rowStride = A.RowStride();

    if( uplo == UPPER )
    {
        T* localBuffer = A.LocalBuffer();
        const int ldim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            int j = rowShift + jLocal*rowStride;
            int lastRow = ( side==LEFT ? j-offset : j-offset+height-width );
            int boundary = std::min( lastRow+1, height );
            int numRows = RawLocalLength( boundary, colShift, colStride );
            T* col = &localBuffer[jLocal*ldim];
            for( int iLocal=0; iLocal<numRows; ++iLocal )
                col[iLocal] *= alpha;
        }
    }
    else
    {
        T* localBuffer = A.LocalBuffer();
        const int ldim = A.LocalLDim();
#ifdef HAVE_OPENMP
        #pragma omp parallel for
#endif
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            int j = rowShift + jLocal*rowStride;
            int firstRow =
                ( side==LEFT ? std::max(j-offset,0)
                             : std::max(j-offset+height-width,0) );
            int numZeroRows = RawLocalLength( firstRow, colShift, colStride );
            T* col = &localBuffer[numZeroRows+jLocal*ldim];
            for( int iLocal=0; iLocal<(localHeight-numZeroRows); ++iLocal )
                col[iLocal] *= alpha;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
