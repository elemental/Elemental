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
