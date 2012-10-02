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
Scale( T alpha, Matrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("Scale");
#endif
    if( alpha != T(1) )
    {
        if( alpha == T(0) )
            for( int j=0; j<X.Width(); ++j )
                for( int i=0; i<X.Height(); ++i )
                    X.Set(i,j,0);
        else
            for( int j=0; j<X.Width(); ++j )
                blas::Scal( X.Height(), alpha, X.Buffer(0,j), 1 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Scale( typename Base<T>::type alpha, Matrix<T>& X )
{ Scale( T(alpha), X ); }

template<typename T>
inline void
Scal( T alpha, Matrix<T>& X )
{ Scale( alpha, X ); }

template<typename T>
inline void
Scal( typename Base<T>::type alpha, Matrix<T>& X )
{ Scale( T(alpha), X ); }

template<typename T,Distribution U,Distribution V>
inline void
Scale( T alpha, DistMatrix<T,U,V>& A )
{ Scale( alpha, A.LocalMatrix() ); }

template<typename T,Distribution U,Distribution V>
inline void
Scale( typename Base<T>::type alpha, DistMatrix<T,U,V>& A )
{ Scale( T(alpha), A.LocalMatrix() ); }

template<typename T,Distribution U,Distribution V>
inline void
Scal( T alpha, DistMatrix<T,U,V>& A )
{ Scale( alpha, A ); }

template<typename T,Distribution U,Distribution V>
inline void
Scal( typename Base<T>::type alpha, DistMatrix<T,U,V>& A )
{ Scale( T(alpha), A ); }

} // namespace elem
