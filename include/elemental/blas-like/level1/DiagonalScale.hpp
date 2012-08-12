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
DiagonalScale
( LeftOrRight side, Orientation orientation,
  const Matrix<T>& d, Matrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("DiagonalScale");
#endif
    const int m = X.Height();
    const int n = X.Width();
    const int ldim = X.LDim();
    if( side == LEFT )
    {
        for( int i=0; i<m; ++i )
        {
            const T delta = d.Get(i,0);
            T* XBuffer = X.Buffer(i,0);
            if( orientation == ADJOINT )
                for( int j=0; j<n; ++j )
                    XBuffer[j*ldim] *= Conj(delta);
            else
                for( int j=0; j<n; ++j )
                    XBuffer[j*ldim] *= delta;
        }
    }
    else
    {
        for( int j=0; j<n; ++j )
        {
            const T delta = d.Get(j,0);
            T* XBuffer = X.Buffer(0,j);
            if( orientation == ADJOINT )
                for( int i=0; i<m; ++i )
                    XBuffer[i] *= Conj(delta);
            else
                for( int i=0; i<m; ++i )
                    XBuffer[i] *= delta;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
DiagonalScale
( LeftOrRight side, Orientation orientation,
  const Matrix<typename Base<T>::type>& d, Matrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("DiagonalScale");
#endif
    typedef typename Base<T>::type R;

    const int m = X.Height();
    const int n = X.Width();
    const int ldim = X.LDim();
    if( side == LEFT )
    {
        for( int i=0; i<m; ++i )
        {
            const R delta = d.Get(i,0);
            T* XBuffer = X.Buffer(i,0);
            for( int j=0; j<n; ++j )
                XBuffer[j*ldim] *= delta;
        }
    }
    else
    {
        for( int j=0; j<n; ++j )
        {
            const R delta = d.Get(j,0);
            T* XBuffer = X.Buffer(0,j);
            for( int i=0; i<m; ++i )
                XBuffer[i] *= delta;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
DiagonalScale
( LeftOrRight side, Orientation orientation,
  const DistMatrix<T,U,V>& d, DistMatrix<T,W,Z>& X )
{
#ifndef RELEASE
    PushCallStack("DiagonalScale");
#endif
    if( side == LEFT )
    {
        if( U == W && V == STAR && d.ColAlignment() == X.ColAlignment() )
        {
            DiagonalScale
            ( LEFT, orientation, d.LockedLocalMatrix(), X.LocalMatrix() );
        }
        else
        {
            DistMatrix<T,W,STAR> d_W_STAR( X.Grid() );
            d_W_STAR = d;
            DiagonalScale
            ( LEFT, orientation,
              d_W_STAR.LockedLocalMatrix(), X.LocalMatrix() );
        }
    }
    else
    {
        if( U == Z && V == STAR && d.ColAlignment() == X.RowAlignment() )
        {
            DiagonalScale
            ( RIGHT, orientation, d.LockedLocalMatrix(), X.LocalMatrix() );
        }
        else
        {
            DistMatrix<T,Z,STAR> d_Z_STAR( X.Grid() );
            d_Z_STAR = d;
            DiagonalScale
            ( RIGHT, orientation,
              d_Z_STAR.LockedLocalMatrix(), X.LocalMatrix() );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
DiagonalScale
( LeftOrRight side, Orientation orientation,
  const DistMatrix<typename Base<T>::type,U,V>& d, DistMatrix<T,W,Z>& X )
{
#ifndef RELEASE
    PushCallStack("DiagonalScale");
#endif
    typedef typename Base<T>::type R;

    if( side == LEFT )
    {
        if( U == W && V == STAR && d.ColAlignment() == X.ColAlignment() )
        {
            DiagonalScale
            ( LEFT, orientation, d.LockedLocalMatrix(), X.LocalMatrix() );
        }
        else
        {
            DistMatrix<R,W,STAR> d_W_STAR( X.Grid() );
            d_W_STAR = d;
            DiagonalScale
            ( LEFT, orientation,
              d_W_STAR.LockedLocalMatrix(), X.LocalMatrix() );
        }
    }
    else
    {
        if( U == Z && V == STAR && d.ColAlignment() == X.RowAlignment() )
        {
            DiagonalScale
            ( RIGHT, orientation, d.LockedLocalMatrix(), X.LocalMatrix() );
        }
        else
        {
            DistMatrix<R,Z,STAR> d_Z_STAR( X.Grid() );
            d_Z_STAR = d;
            DiagonalScale
            ( RIGHT, orientation,
              d_Z_STAR.LockedLocalMatrix(), X.LocalMatrix() );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
