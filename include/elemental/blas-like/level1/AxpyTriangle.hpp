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
AxpyTriangle( UpperOrLower uplo, T alpha, const Matrix<T>& X, Matrix<T>& Y )
{
#ifndef RELEASE
    PushCallStack("AxpyTriangle");
    if( X.Height() != X.Width() || Y.Height() != Y.Width() || 
        X.Height() != Y.Height() )
        throw std::logic_error("Nonconformal AxpyTriangle");
#endif
    if( uplo == UPPER )
    {
        for( int j=0; j<X.Width(); ++j )
            blas::Axpy( j+1, alpha, X.LockedBuffer(0,j), 1, Y.Buffer(0,j), 1 );
    }
    else
    {
        const int n = X.Height();
        for( int j=0; j<X.Width(); ++j )
            blas::Axpy( n-j, alpha, X.LockedBuffer(j,j), 1, Y.Buffer(j,j), 1 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
AxpyTriangle
( UpperOrLower uplo, typename Base<T>::type alpha, 
  const Matrix<T>& X, Matrix<T>& Y )
{ AxpyTriangle( uplo, T(alpha), X, Y ); }

template<typename T,Distribution U,Distribution V>
inline void
AxpyTriangle
( UpperOrLower uplo, T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )
{
#ifndef RELEASE
    PushCallStack("AxpyTriangle");
    if( X.Grid() != Y.Grid() )
        throw std::logic_error
        ("X and Y must be distributed over the same grid");
    if( X.Height() != X.Width() || Y.Height() != Y.Width() || 
        X.Height() != Y.Height() )
        throw std::logic_error("Nonconformal AxpyTriangle");
#endif
    if( X.ColAlignment() == Y.ColAlignment() &&
        X.RowAlignment() == Y.RowAlignment() )
    {
        const int localHeight = X.LocalHeight();
        const int localWidth = X.LocalWidth();
        const int colShift = X.ColShift();
        const int rowShift = X.RowShift();
        const int colStride = X.ColStride();
        const int rowStride = X.RowStride();
        const T* XBuffer = X.LockedLocalBuffer();
        T* YBuffer = Y.LocalBuffer();
        const int XLDim = X.LocalLDim();
        const int YLDim = Y.LocalLDim();
        if( uplo == UPPER )
        {
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const int j = rowShift + jLocal*rowStride;        
                const int localHeightAbove = 
                    LocalLength( j+1, colShift, colStride );
                blas::Axpy
                ( localHeightAbove, alpha, 
                  &XBuffer[jLocal*XLDim], 1, &YBuffer[jLocal*YLDim], 1 );
            }
        }
        else
        {
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const int j = rowShift + jLocal*rowStride;
                const int localHeightAbove = 
                    LocalLength( j, colShift, colStride );
                const int localHeightBelow = localHeight - localHeightAbove;
                blas::Axpy
                ( localHeightBelow, alpha, 
                  &XBuffer[localHeightAbove+jLocal*XLDim], 1,
                  &YBuffer[localHeightAbove+jLocal*YLDim], 1 );
            }
        }
    }
    else
    {
        DistMatrix<T,U,V> XCopy( X.Grid() );
        XCopy.AlignWith( Y );
        XCopy = X;
        AxpyTriangle( uplo, alpha, XCopy.LockedLocalMatrix(), Y.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
AxpyTriangle
( UpperOrLower uplo, typename Base<T>::type alpha,
  const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )
{ AxpyTriangle( uplo, T(alpha), X, Y ); }

} // namespace elem
