/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
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
