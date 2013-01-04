/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

template<typename R>
inline void
DiscreteFourier( int n, Matrix<Complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("DiscreteFourier");
#endif
    A.ResizeTo( n, n );
    MakeDiscreteFourier( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R,Distribution U,Distribution V>
inline void
DiscreteFourier( int n, DistMatrix<Complex<R>,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("DiscreteFourier");
#endif
    A.ResizeTo( n, n );
    MakeDiscreteFourier( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
MakeDiscreteFourier( Matrix<Complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("MakeDiscreteFourier");
#endif
    typedef Complex<R> F;

    const int m = A.Height();
    const int n = A.Width();
    if( m != n )
        throw std::logic_error("Cannot make a non-square DFT matrix");

    const R pi = 4*Atan( R(1) );
    const F nSqrt = Sqrt( R(n) );
    for( int j=0; j<n; ++j )
    {
        for( int i=0; i<m; ++i )
        {
            const R theta = -2*pi*i*j/n;
            A.Set( i, j, Complex<R>(Cos(theta),Sin(theta))/nSqrt );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R,Distribution U,Distribution V>
inline void
MakeDiscreteFourier( DistMatrix<Complex<R>,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("MakeDiscreteFourier");
#endif
    typedef Complex<R> F;

    const int m = A.Height();
    const int n = A.Width();
    if( m != n )
        throw std::logic_error("Cannot make a non-square DFT matrix");

    const R pi = 4*Atan( R(1) );
    const F nSqrt = Sqrt( R(n) );
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();
    const int colStride = A.ColStride();
    const int rowStride = A.RowStride();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const int j = rowShift + jLocal*rowStride;
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = colShift + iLocal*colStride;
            A.SetLocal( iLocal, jLocal, Exp(-2*pi*i*j/n)/nSqrt );

            const R theta = -2*pi*i*j/n;
            const Complex<R> alpha( Cos(theta), Sin(theta) );
            A.SetLocal( iLocal, jLocal, alpha/nSqrt );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
