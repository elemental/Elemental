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
