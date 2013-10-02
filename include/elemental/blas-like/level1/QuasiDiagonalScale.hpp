/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_QUASIDIAGONALSCALE_HPP
#define ELEM_BLAS_QUASIDIAGONALSCALE_HPP

#include "elemental/blas-like/level1/Symmetric2x2Scale.hpp"

namespace elem {

template<typename F,typename FMain>
inline void
QuasiDiagonalScale
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const Matrix<FMain>& d, const Matrix<F>& dSub, 
  Matrix<F>& X, bool conjugated=false )
{
#ifndef RELEASE
    CallStackEntry entry("QuasiDiagonalSolve");
#endif
    const Int m = X.Height();
    const Int n = X.Width();
    Matrix<F> D( 2, 2 );
    if( side == LEFT && uplo == LOWER )
    {
        Int i=0;
        while( i < m )
        {
            Int nb;
            if( i < m-1 && Abs(dSub.Get(i,0)) > 0 )
                nb = 2;
            else
                nb = 1;

            if( nb == 1 )
            {
                auto xRow = View( X, i, 0, nb, n );
                Scale( d.Get(i,0), xRow );
            }
            else
            {
                D.Set(0,0,d.Get(i,0));    
                D.Set(1,1,d.Get(i+1,0));
                D.Set(1,0,dSub.Get(i,0));
                auto XRow = View( X, i, 0, nb, n );
                Symmetric2x2Scale( LEFT, LOWER, D, XRow, conjugated );
            }

            i += nb;
        }
    }
    else if( side == RIGHT && uplo == LOWER )
    {
        Int j=0;
        while( j < n )
        {
            Int nb;
            if( j < n-1 && Abs(dSub.Get(j,0)) > 0 )
                nb = 2;
            else
                nb = 1;

            if( nb == 1 )
            {
                auto xCol = View( X, 0, j, m, nb );
                Scale( d.Get(j,0), xCol );
            }
            else
            {
                D.Set(0,0,d.Get(j,0));    
                D.Set(1,1,d.Get(j+1,0));
                D.Set(1,0,dSub.Get(j,0));
                auto XCol = View( X, 0, j, m, nb );
                Symmetric2x2Scale( RIGHT, LOWER, D, XCol, conjugated );
            }

            j += nb;
        }
    }
    else
        LogicError("This option not yet supported");
}

template<typename F,typename FMain>
inline void
LeftQuasiDiagonalScale
( UpperOrLower uplo, Orientation orientation, 
  const DistMatrix<FMain,MC,STAR> d,
  const DistMatrix<FMain,MC,STAR> dPrev,
  const DistMatrix<FMain,MC,STAR> dNext,
  const DistMatrix<FMain,MC,STAR> dSub,
  const DistMatrix<FMain,MC,STAR> dSubPrev,
  const DistMatrix<FMain,MC,STAR> dSubNext,
        DistMatrix<F>& X,
  const DistMatrix<F>& XPrev,
  const DistMatrix<F>& XNext,
  bool conjugated=false )
{
#ifndef RELEASE
    CallStackEntry entry("LeftQuasiDiagonalScale");
#endif
    if( uplo == UPPER || orientation != NORMAL )
        LogicError("This option not yet supported");
    const Int mLocal = X.LocalHeight();
    const Int colShift = X.ColShift();
    const Int colStride = X.ColStride();
    const Int colAlignPrev = (X.ColAlign()+colStride-1) % colStride;
    const Int colAlignNext = (X.ColAlign()+1) % colStride;
#ifndef RELEASE
    if( d.ColAlign() != X.ColAlign() || dSub.ColAlign() != X.ColAlign() )
        LogicError("data is not properly aligned");
    if( XPrev.ColAlign() != colAlignPrev ||
        dPrev.ColAlign() != colAlignPrev || 
        dSubPrev.ColAlign() != colAlignPrev )
        LogicError("'previous' data is not properly aligned");
    if( XNext.ColAlign() != colAlignNext || 
        dNext.ColAlign() != colAlignNext || 
        dSubNext.ColAlign() != colAlignNext )
        LogicError("'next' data is not properly aligned");
#endif
    const Int colShiftPrev = XPrev.ColShift();
    const Int colShiftNext = XNext.ColShift();
    const Int prevOff = ( colShiftPrev==colShift-1 ? 0 : -1 );
    const Int nextOff = ( colShiftNext==colShift+1 ? 0 : +1 );
    const Int m = X.Height();
    const Int nLocal = X.LocalWidth();

    // It is best to separate the case where colStride is 1
    if( colStride == 1 )
    {
        QuasiDiagonalScale
        ( LEFT, uplo, orientation, d.LockedMatrix(), dSub.LockedMatrix(),
          X.Matrix(), conjugated );
        return;
    }

    Matrix<F> D11( 2, 2 );
    for( Int iLoc=0; iLoc<mLocal; ++iLoc )
    {
        const Int i = colShift + iLoc*colStride;
        const Int iLocPrev = iLoc + prevOff;
        const Int iLocNext = iLoc + nextOff;

        auto x1Loc = View( X.Matrix(), iLoc, 0, 1, nLocal );

        if( i<m-1 && dSub.GetLocal(iLoc,0) != F(0) )
        {
            // Handle 2x2 starting at i
            D11.Set( 0, 0, d.GetLocal(iLoc,0) ); 
            D11.Set( 1, 1, dNext.GetLocal(iLocNext,0) );
            D11.Set( 1, 0, dSub.GetLocal(iLoc,0) );

            auto x1NextLoc = 
                LockedView( XNext.LockedMatrix(), iLocNext, 0, 1, nLocal );
            FirstHalfOfSymmetric2x2Scale
            ( LEFT, LOWER, D11, x1Loc, x1NextLoc, conjugated );
        }
        else if( i>0 && dSubPrev.GetLocal(iLocPrev,0) != F(0) )
        {
            // Handle 2x2 starting at i-1
            D11.Set( 0, 0, dPrev.GetLocal(iLocPrev,0) );
            D11.Set( 1, 1, d.GetLocal(iLoc,0) );
            D11.Set( 1, 0, dSubPrev.GetLocal(iLocPrev,0) );

            auto x1PrevLoc = 
                LockedView( XPrev.LockedMatrix(), iLocPrev, 0, 1, nLocal );
            SecondHalfOfSymmetric2x2Scale
            ( LEFT, LOWER, D11, x1PrevLoc, x1Loc, conjugated );
        }
        else
        {
            // Handle 1x1
            Scale( d.GetLocal(iLoc,0), x1Loc );
        }
    }
}

template<typename F,typename FMain,Distribution U,Distribution V>
inline void
QuasiDiagonalScale
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const DistMatrix<FMain,U,V>& d, const DistMatrix<F,U,V>& dSub, 
  DistMatrix<F>& X, bool conjugated=false )
{
#ifndef RELEASE
    CallStackEntry entry("QuasiDiagonalScale");
#endif
    const Grid& g = X.Grid();
    const Int colAlign = X.ColAlign();
    const Int rowAlign = X.RowAlign();
    const Int colStride = X.ColStride();
    if( side == LEFT )
    {
        DistMatrix<FMain,MC,STAR> d_MC_STAR(g);
        DistMatrix<F,MC,STAR> dSub_MC_STAR(g);
        d_MC_STAR.AlignWith( X );
        dSub_MC_STAR.AlignWith( X );
        d_MC_STAR = d;
        dSub_MC_STAR = dSub;
        if( colStride == 1 )
        {
            QuasiDiagonalScale
            ( side, uplo, orientation, 
              d_MC_STAR.LockedMatrix(), dSub_MC_STAR.LockedMatrix(),
              X.Matrix(), conjugated );
            return;
        }

        DistMatrix<FMain,MC,STAR> dPrev_MC_STAR(g), dNext_MC_STAR(g);
        DistMatrix<F,MC,STAR> dSubPrev_MC_STAR(g), dSubNext_MC_STAR(g);
        DistMatrix<F> XPrev(g), XNext(g);
        const Int colAlignPrev = (colAlign+colStride-1) % colStride;
        const Int colAlignNext = (colAlign+1) % colStride;
        dPrev_MC_STAR.AlignCols( colAlignPrev );
        dNext_MC_STAR.AlignCols( colAlignNext );
        dSubPrev_MC_STAR.AlignCols( colAlignPrev );
        dSubNext_MC_STAR.AlignCols( colAlignNext );
        XPrev.Align( colAlignPrev, rowAlign );
        XNext.Align( colAlignNext, rowAlign );
        dPrev_MC_STAR = d;
        dNext_MC_STAR = d;
        dSubPrev_MC_STAR = dSub;
        dSubNext_MC_STAR = dSub;
        XPrev = X;
        XNext = X;
        LeftQuasiDiagonalScale
        ( uplo, orientation, 
          d_MC_STAR, dPrev_MC_STAR, dNext_MC_STAR,
          dSub_MC_STAR, dSubPrev_MC_STAR, dSubNext_MC_STAR,
          X, XPrev, XNext, conjugated );
    }
    else
        LogicError("Not yet written");
}

} // namespace elem

#endif // ifndef ELEM_BLAS_QUASIDIAGONALSCALE_HPP
