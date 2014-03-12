/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_DIAGONALSCALE_HPP
#define ELEM_DIAGONALSCALE_HPP

namespace elem {

template<typename T,typename TDiag>
inline void
DiagonalScale
( LeftOrRight side, Orientation orientation,
  const Matrix<TDiag>& d, Matrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalScale"))
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

template<typename T,typename TDiag,
         Dist U,Dist V,
         Dist W,Dist Z>
inline void
DiagonalScale
( LeftOrRight side, Orientation orientation,
  const DistMatrix<TDiag,U,V>& d, DistMatrix<T,W,Z>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalScale"))
    if( side == LEFT )
    {
        DistMatrix<TDiag,W,STAR> d_W_STAR( X.Grid() );
        if( U == W && V == STAR && d.ColAlign() == X.ColAlign() )
        {
            d_W_STAR = LockedView( d );
        }
        else
        {
            d_W_STAR.AlignWith( X );
            d_W_STAR = d;
        }
        DiagonalScale( LEFT, orientation, d_W_STAR.LockedMatrix(), X.Matrix() );
    }
    else
    {
        DistMatrix<TDiag,Z,STAR> d_Z_STAR( X.Grid() );
        if( U == Z && V == STAR && d.ColAlign() == X.RowAlign() )
        {
            d_Z_STAR = LockedView( d );
        }
        else
        {
            d_Z_STAR.AlignWith( X );
            d_Z_STAR = d;
        }
        DiagonalScale
        ( RIGHT, orientation, d_Z_STAR.LockedMatrix(), X.Matrix() );
    }
}

} // namespace elem

#endif // ifndef ELEM_DIAGONALSCALE_HPP
