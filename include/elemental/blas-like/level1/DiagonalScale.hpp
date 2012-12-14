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
