/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_DIAGONALSCALE_HPP
#define ELEM_BLAS_DIAGONALSCALE_HPP

namespace elem {

template<typename T>
inline void
DiagonalScale
( LeftOrRight side, Orientation orientation,
  const Matrix<T>& d, Matrix<T>& X )
{
#ifndef RELEASE
    CallStackEntry entry("DiagonalScale");
#endif
    const Int m = X.Height();
    const Int n = X.Width();
    const Int ldim = X.LDim();
    if( side == LEFT )
    {
        for( Int i=0; i<m; ++i )
        {
            const T delta = d.Get(i,0);
            T* XBuffer = X.Buffer(i,0);
            if( orientation == ADJOINT )
                for( Int j=0; j<n; ++j )
                    XBuffer[j*ldim] *= Conj(delta);
            else
                for( Int j=0; j<n; ++j )
                    XBuffer[j*ldim] *= delta;
        }
    }
    else
    {
        for( Int j=0; j<n; ++j )
        {
            const T delta = d.Get(j,0);
            T* XBuffer = X.Buffer(0,j);
            if( orientation == ADJOINT )
                for( Int i=0; i<m; ++i )
                    XBuffer[i] *= Conj(delta);
            else
                for( Int i=0; i<m; ++i )
                    XBuffer[i] *= delta;
        }
    }
}

template<typename T>
inline void
DiagonalScale
( LeftOrRight side, Orientation orientation,
  const Matrix<BASE(T)>& d, Matrix<T>& X )
{
#ifndef RELEASE
    CallStackEntry entry("DiagonalScale");
#endif
    typedef BASE(T) R;

    const Int m = X.Height();
    const Int n = X.Width();
    const Int ldim = X.LDim();
    if( side == LEFT )
    {
        for( Int i=0; i<m; ++i )
        {
            const R delta = d.Get(i,0);
            T* XBuffer = X.Buffer(i,0);
            for( Int j=0; j<n; ++j )
                XBuffer[j*ldim] *= delta;
        }
    }
    else
    {
        for( Int j=0; j<n; ++j )
        {
            const R delta = d.Get(j,0);
            T* XBuffer = X.Buffer(0,j);
            for( Int i=0; i<m; ++i )
                XBuffer[i] *= delta;
        }
    }
}

template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
DiagonalScale
( LeftOrRight side, Orientation orientation,
  const DistMatrix<T,U,V>& d, DistMatrix<T,W,Z>& X )
{
#ifndef RELEASE
    CallStackEntry entry("DiagonalScale");
#endif
    if( side == LEFT )
    {
        if( U == W && V == STAR && d.ColAlign() == X.ColAlign() )
        {
            DiagonalScale( LEFT, orientation, d.LockedMatrix(), X.Matrix() );
        }
        else
        {
            DistMatrix<T,W,STAR> d_W_STAR( X.Grid() );
            d_W_STAR = d;
            DiagonalScale
            ( LEFT, orientation,
              d_W_STAR.LockedMatrix(), X.Matrix() );
        }
    }
    else
    {
        if( U == Z && V == STAR && d.ColAlign() == X.RowAlign() )
        {
            DiagonalScale
            ( RIGHT, orientation, d.LockedMatrix(), X.Matrix() );
        }
        else
        {
            DistMatrix<T,Z,STAR> d_Z_STAR( X.Grid() );
            d_Z_STAR = d;
            DiagonalScale
            ( RIGHT, orientation, d_Z_STAR.LockedMatrix(), X.Matrix() );
        }
    }
}

template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
DiagonalScale
( LeftOrRight side, Orientation orientation,
  const DistMatrix<BASE(T),U,V>& d, DistMatrix<T,W,Z>& X )
{
#ifndef RELEASE
    CallStackEntry entry("DiagonalScale");
#endif
    typedef BASE(T) R;

    if( side == LEFT )
    {
        if( U == W && V == STAR && d.ColAlign() == X.ColAlign() )
        {
            DiagonalScale( LEFT, orientation, d.LockedMatrix(), X.Matrix() );
        }
        else
        {
            DistMatrix<R,W,STAR> d_W_STAR( X.Grid() );
            d_W_STAR = d;
            DiagonalScale
            ( LEFT, orientation, d_W_STAR.LockedMatrix(), X.Matrix() );
        }
    }
    else
    {
        if( U == Z && V == STAR && d.ColAlign() == X.RowAlign() )
        {
            DiagonalScale( RIGHT, orientation, d.LockedMatrix(), X.Matrix() );
        }
        else
        {
            DistMatrix<R,Z,STAR> d_Z_STAR( X.Grid() );
            d_Z_STAR = d;
            DiagonalScale
            ( RIGHT, orientation, d_Z_STAR.LockedMatrix(), X.Matrix() );
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_BLAS_DIAGONALSCALE_HPP
