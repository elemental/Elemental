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

template<typename F> 
inline void
Legendre( int n, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Legendre");
#endif
    A.ResizeTo( n, n );
    MakeLegendre( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F,Distribution U,Distribution V> 
inline void
Legendre( int n, DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("Legendre");
#endif
    A.ResizeTo( n, n );
    MakeLegendre( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
MakeLegendre( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("MakeLegendre");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make a non-square matrix Legendre");
    MakeZeros( A );

    const int n = A.Width();
    for( int j=0; j<n-1; ++j )
    {
        const F gamma = F(1) / Pow( F(2)*(j+1), F(2) );
        const F beta = F(1) / (2*Sqrt(F(1)-gamma));
        A.Set( j+1, j, beta );
        A.Set( j, j+1, beta );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F,Distribution U,Distribution V>
inline void
MakeLegendre( DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("MakeLegendre");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot make a non-square matrix Legendre");
    MakeZeros( A );

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
            if( j == i+1 || j == i-1 )
            {
                const int k = std::max( i, j );
                const F gamma = F(1) / Pow( F(2)*k, F(2) );
                const F beta = F(1) / (2*Sqrt(F(1)-gamma));
                A.SetLocal( iLocal, jLocal, beta );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
