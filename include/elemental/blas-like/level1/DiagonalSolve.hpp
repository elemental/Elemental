/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_DIAGONALSOLVE_HPP
#define ELEM_DIAGONALSOLVE_HPP

namespace elem {

template<typename F,typename FDiag>
inline void
DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const Matrix<FDiag>& d, Matrix<F>& X, bool checkIfSingular=true )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalSolve"))
    const Int m = X.Height();
    const Int n = X.Width();
    const Int ldim = X.LDim();
    if( side == LEFT )
    {
        for( Int i=0; i<m; ++i )
        {
            const F delta = d.Get(i,0);
            if( checkIfSingular && delta == F(0) )
                throw SingularMatrixException();
            const F deltaInv = F(1)/delta;
            F* XBuffer = X.Buffer(i,0);
            if( orientation == ADJOINT )
                for( Int j=0; j<n; ++j )
                    XBuffer[j*ldim] *= Conj(deltaInv);
            else
                for( Int j=0; j<n; ++j )
                    XBuffer[j*ldim] *= deltaInv;
        }
    }
    else
    {
        for( Int j=0; j<n; ++j )
        {
            const F delta = d.Get(j,0);
            if( checkIfSingular && delta == F(0) )
                throw SingularMatrixException();
            const F deltaInv = F(1)/delta;
            F* XBuffer = X.Buffer(0,j);
            if( orientation == ADJOINT )
                for( Int i=0; i<m; ++i )
                    XBuffer[i] *= Conj(deltaInv);
            else
                for( Int i=0; i<m; ++i )
                    XBuffer[i] *= deltaInv;
        }
    }
}

template<typename F,typename FDiag,
         Dist U,Dist V,
         Dist W,Dist Z>
inline void
DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const DistMatrix<FDiag,U,V>& d, DistMatrix<F,W,Z>& X,
  bool checkIfSingular=true )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalSolve"))
    if( side == LEFT )
    {
        if( U == W && V == STAR && d.ColAlign() == X.ColAlign() )
        {
            DiagonalSolve
            ( LEFT, orientation, d.LockedMatrix(), X.Matrix(),
              checkIfSingular );
        }
        else
        {
            DistMatrix<FDiag,W,STAR> d_W_STAR( X.Grid() );
            d_W_STAR = d;
            DiagonalSolve
            ( LEFT, orientation,
              d_W_STAR.LockedMatrix(), X.Matrix(), checkIfSingular );
        }
    }
    else
    {
        if( U == Z && V == STAR && d.ColAlign() == X.RowAlign() )
        {
            DiagonalSolve
            ( RIGHT, orientation, d.LockedMatrix(), X.Matrix(),
              checkIfSingular );
        }
        else
        {
            DistMatrix<FDiag,Z,STAR> d_Z_STAR( X.Grid() );
            d_Z_STAR = d;
            DiagonalSolve
            ( RIGHT, orientation,
              d_Z_STAR.LockedMatrix(), X.Matrix(), checkIfSingular );
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_DIAGONALSOLVE_HPP
