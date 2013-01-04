/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

template<typename T>
inline void
MakeTrapezoidal
( LeftOrRight side, UpperOrLower uplo, int offset, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("MakeTrapezoidal");
#endif
    const int height = A.Height();
    const int width = A.Width();
    const int ldim = A.LDim();
    T* buffer = A.Buffer();

    if( uplo == LOWER )
    {
        if( side == LEFT )
        {
#ifdef HAVE_OPENMP
            #pragma omp parallel for
#endif
            for( int j=std::max(0,offset+1); j<width; ++j )
            {
                const int lastZeroRow = j-offset-1;
                const int numZeroRows = std::min( lastZeroRow+1, height );
                MemZero( &buffer[j*ldim], numZeroRows );
            }
        }
        else
        {
#ifdef HAVE_OPENMP
            #pragma omp parallel for
#endif
            for( int j=std::max(0,offset-height+width+1); j<width; ++j )
            {
                const int lastZeroRow = j-offset+height-width-1;
                const int numZeroRows = std::min( lastZeroRow+1, height );
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
            for( int j=0; j<width; ++j )
            {
                const int firstZeroRow = std::max(j-offset+1,0);
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
            for( int j=0; j<width; ++j )
            {
                const int firstZeroRow = std::max(j-offset+height-width+1,0);
                if( firstZeroRow < height )
                    MemZero
                    ( &buffer[firstZeroRow+j*ldim], height-firstZeroRow );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
MakeTrapezoidal
( LeftOrRight side, UpperOrLower uplo, int offset,
  DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("MakeTrapezoidal");
#endif
    const int height = A.Height();
    const int width = A.Width();
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
            const int lastZeroRow =
                ( side==LEFT ? j-offset-1
                             : j-offset+height-width-1 );
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
            const int firstZeroRow =
                ( side==LEFT ? std::max(j-offset+1,0)
                             : std::max(j-offset+height-width+1,0) );
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
