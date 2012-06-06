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

    const R pi = 4*Atan( (R)1 );
    const F nSqrt = Sqrt( (R)n );
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

    const R pi = 4*Atan( (R)1 );
    const F nSqrt = Sqrt( (R)n );
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
