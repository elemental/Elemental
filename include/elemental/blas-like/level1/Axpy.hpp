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
Axpy( T alpha, const Matrix<T>& X, Matrix<T>& Y )
{
#ifndef RELEASE
    PushCallStack("Axpy");
#endif
    // If X and Y are vectors, we can allow one to be a column and the other
    // to be a row. Otherwise we force X and Y to be the same dimension.
    if( (X.Height()==1 || X.Width()==1) && (Y.Height()==1 || Y.Width()==1) )
    {
        const unsigned XLength = ( X.Width()==1 ? X.Height() : X.Width() );
#ifndef RELEASE
        const unsigned YLength = ( Y.Width()==1 ? Y.Height() : Y.Width() );
        if( XLength != YLength )
            throw std::logic_error("Nonconformal Axpy");
#endif
        if( X.Width()==1 && Y.Width()==1 )
        {
            blas::Axpy
            ( XLength, alpha, X.LockedBuffer(0,0), 1, Y.Buffer(0,0), 1 );
        }
        else if( X.Width()==1 )
        {
            blas::Axpy
            ( XLength, alpha, X.LockedBuffer(0,0), 1, Y.Buffer(0,0), Y.LDim() );
        }
        else if( Y.Width()==1 )
        {
            blas::Axpy
            ( XLength, alpha, X.LockedBuffer(0,0), X.LDim(), Y.Buffer(0,0), 1 );
        }
        else
        {
            blas::Axpy
            ( XLength, alpha,
              X.LockedBuffer(0,0), X.LDim(), Y.Buffer(0,0), Y.LDim() );
        }
    }
    else
    {
#ifndef RELEASE
        if( X.Height() != Y.Height() || X.Width() != Y.Width() )
            throw std::logic_error("Nonconformal Axpy");
#endif
        if( X.Width() <= X.Height() )
        {
            for( int j=0; j<X.Width(); ++j )
            {
                blas::Axpy
                ( X.Height(), alpha, X.LockedBuffer(0,j), 1, Y.Buffer(0,j), 1 );
            }
        }
        else
        {
            for( int i=0; i<X.Height(); ++i )
            {
                blas::Axpy
                ( X.Width(), alpha, X.LockedBuffer(i,0), X.LDim(),
                                    Y.Buffer(i,0),       Y.LDim() );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Axpy( typename Base<T>::type alpha, const Matrix<T>& X, Matrix<T>& Y )
{ Axpy( T(alpha), X, Y ); }

template<typename T,Distribution U,Distribution V>
inline void
Axpy( T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )
{
#ifndef RELEASE
    PushCallStack("Axpy");
    if( X.Grid() != Y.Grid() )
        throw std::logic_error
        ("X and Y must be distributed over the same grid");
#endif
    if( X.ColAlignment() == Y.ColAlignment() &&
        X.RowAlignment() == Y.RowAlignment() )
    {
        Axpy( alpha, X.LockedLocalMatrix(), Y.LocalMatrix() );
    }
    else
    {
        DistMatrix<T,U,V> XCopy( X.Grid() );
        XCopy.AlignWith( Y );
        XCopy = X;
        Axpy( alpha, XCopy.LockedLocalMatrix(), Y.LocalMatrix() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
Axpy
( typename Base<T>::type alpha,
  const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )
{ Axpy( T(alpha), X, Y ); }

} // namespace elem
