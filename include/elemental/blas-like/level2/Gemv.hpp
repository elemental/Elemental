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

#include "./Gemv/N.hpp"
#include "./Gemv/T.hpp"

namespace elem {

template<typename T>
inline void
Gemv
( Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Gemv");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
    {
        std::ostringstream msg;
        msg << "x and y must be vectors:\n"
            << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
            << "  y ~ " << y.Height() << " x " << y.Width();
        throw std::logic_error( msg.str().c_str() );
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( orientation == NORMAL )
    {
        if( A.Height() != yLength || A.Width() != xLength )
        {
            std::ostringstream msg;
            msg << "A must conform with x and y:\n"
                << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
                << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
                << "  y ~ " << y.Height() << " x " << y.Width();
            throw std::logic_error( msg.str().c_str() );
        }
    }
    else
    {
        if( A.Width() != yLength || A.Height() != xLength )
        {
            std::ostringstream msg;
            msg << "A must conform with x and y:\n"
                << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
                << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
                << "  y ~ " << y.Height() << " x " << y.Width();
            throw std::logic_error( msg.str().c_str() );
        }
    }
#endif
    const char transChar = OrientationToChar( orientation );
    const int m = A.Height();
    const int n = A.Width();
    const int k = ( transChar == 'N' ? n : m );
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    if( k != 0 )
    {
        blas::Gemv
        ( transChar, m, n,
          alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx,
          beta,  y.Buffer(), incy );
    }
    else
    {
        Scale( beta, y );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Gemv
( Orientation orientation,
  T alpha, const DistMatrix<T>& A, 
           const DistMatrix<T>& x,
  T beta,        DistMatrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Gemv");
#endif
    if( orientation == NORMAL )
        internal::GemvN( alpha, A, x, beta, y );
    else
        internal::GemvT( orientation, alpha, A, x, beta, y );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Gemv
( Orientation orientation,
  T alpha, const DistMatrix<T>& A, 
           const DistMatrix<T,VC,STAR>& x,
  T beta,        DistMatrix<T,VC,STAR>& y )
{
#ifndef RELEASE
    PushCallStack("Gemv");
#endif
    if( orientation == NORMAL )
        internal::GemvN( alpha, A, x, beta, y );
    else
        internal::GemvT( orientation, alpha, A, x, beta, y );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

