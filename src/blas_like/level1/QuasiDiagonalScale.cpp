/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F,typename FMain>
void QuasiDiagonalScale
( LeftOrRight side, UpperOrLower uplo, 
  const Matrix<FMain>& d, const Matrix<F>& dSub, 
  Matrix<F>& X, bool conjugated )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiDiagonalScale"))
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

template<typename F,typename FMain,Dist U,Dist V>
void LeftQuasiDiagonalScale
( UpperOrLower uplo,
  const DistMatrix<FMain,U,STAR>& d,
  const DistMatrix<FMain,U,STAR>& dPrev,
  const DistMatrix<FMain,U,STAR>& dNext,
  const DistMatrix<F,    U,STAR>& dSub,
  const DistMatrix<F,    U,STAR>& dSubPrev,
  const DistMatrix<F,    U,STAR>& dSubNext,
        DistMatrix<F,U,V>& X,
  const DistMatrix<F,U,V>& XPrev,
  const DistMatrix<F,U,V>& XNext,
  bool conjugated )
{
    DEBUG_ONLY(CallStackEntry cse("LeftQuasiDiagonalScale"))
    if( uplo == UPPER )
        LogicError("This option not yet supported");
    const Int m = X.Height();
    const Int mLocal = X.LocalHeight();
    const Int nLocal = X.LocalWidth();
    const Int colStride = X.ColStride();
    DEBUG_ONLY(
        const Int colAlignPrev = Mod(X.ColAlign()+1,colStride);
        const Int colAlignNext = Mod(X.ColAlign()-1,colStride);
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
    )
    const Int prevOff = ( XPrev.ColShift()==X.ColShift()-1 ? 0 : -1 );
    const Int nextOff = ( XNext.ColShift()==X.ColShift()+1 ? 0 : +1 );
    if( !X.Participating() )
        return;

    // It is best to separate the case where colStride is 1
    if( colStride == 1 )
    {
        QuasiDiagonalScale
        ( LEFT, uplo, d.LockedMatrix(), dSub.LockedMatrix(), X.Matrix(),
          conjugated );
        return;
    }

    Matrix<F> D11( 2, 2 );
    for( Int iLoc=0; iLoc<mLocal; ++iLoc )
    {
        const Int i = X.GlobalRow(iLoc);
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

template<typename F,typename FMain,Dist U,Dist V>
void RightQuasiDiagonalScale
( UpperOrLower uplo,
  const DistMatrix<FMain,V,STAR>& d,
  const DistMatrix<FMain,V,STAR>& dPrev,
  const DistMatrix<FMain,V,STAR>& dNext,
  const DistMatrix<F,    V,STAR>& dSub,
  const DistMatrix<F,    V,STAR>& dSubPrev,
  const DistMatrix<F,    V,STAR>& dSubNext,
        DistMatrix<F,U,V>& X,
  const DistMatrix<F,U,V>& XPrev,
  const DistMatrix<F,U,V>& XNext,
  bool conjugated )
{
    DEBUG_ONLY(CallStackEntry cse("LeftQuasiDiagonalScale"))
    if( uplo == UPPER )
        LogicError("This option not yet supported");
    const Int n = X.Width();
    const Int mLocal = X.LocalHeight();
    const Int nLocal = X.LocalWidth();
    const Int rowStride = X.RowStride();
    DEBUG_ONLY(
        const Int rowAlignPrev = Mod(X.RowAlign()+1,rowStride);
        const Int rowAlignNext = Mod(X.RowAlign()-1,rowStride);
        if( d.ColAlign() != X.RowAlign() || dSub.ColAlign() != X.RowAlign() )
            LogicError("data is not properly aligned");
        if( XPrev.RowAlign() != rowAlignPrev ||
            dPrev.ColAlign() != rowAlignPrev || 
            dSubPrev.ColAlign() != rowAlignPrev )
            LogicError("'previous' data is not properly aligned");
        if( XNext.RowAlign() != rowAlignNext || 
            dNext.ColAlign() != rowAlignNext || 
            dSubNext.ColAlign() != rowAlignNext )
            LogicError("'next' data is not properly aligned");
    )
    const Int prevOff = ( XPrev.RowShift()==X.RowShift()-1 ? 0 : -1 );
    const Int nextOff = ( XNext.RowShift()==X.RowShift()+1 ? 0 : +1 );
    if( !X.Participating() )
        return;

    // It is best to separate the case where rowStride is 1
    if( rowStride == 1 )
    {
        QuasiDiagonalScale
        ( LEFT, uplo, d.LockedMatrix(), dSub.LockedMatrix(), X.Matrix(),
          conjugated );
        return;
    }

    Matrix<F> D11( 2, 2 );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        const Int j = X.GlobalCol(jLoc);
        const Int jLocPrev = jLoc + prevOff;
        const Int jLocNext = jLoc + nextOff;

        auto x1Loc = View( X.Matrix(), 0, jLoc, mLocal, 1 );

        if( j<n-1 && dSub.GetLocal(jLoc,0) != F(0) )
        {
            // Handle 2x2 starting at j
            D11.Set( 0, 0, d.GetLocal(jLoc,0) ); 
            D11.Set( 1, 1, dNext.GetLocal(jLocNext,0) );
            D11.Set( 1, 0, dSub.GetLocal(jLoc,0) );

            auto x1NextLoc = 
                LockedView( XNext.LockedMatrix(), 0, jLocNext, mLocal, 1 );
            FirstHalfOfSymmetric2x2Scale
            ( RIGHT, LOWER, D11, x1Loc, x1NextLoc, conjugated );
        }
        else if( j>0 && dSubPrev.GetLocal(jLocPrev,0) != F(0) )
        {
            // Handle 2x2 starting at j-1
            D11.Set( 0, 0, dPrev.GetLocal(jLocPrev,0) );
            D11.Set( 1, 1, d.GetLocal(jLoc,0) );
            D11.Set( 1, 0, dSubPrev.GetLocal(jLocPrev,0) );

            auto x1PrevLoc = 
                LockedView( XPrev.LockedMatrix(), 0, jLocPrev, mLocal, 1 );
            SecondHalfOfSymmetric2x2Scale
            ( RIGHT, LOWER, D11, x1PrevLoc, x1Loc, conjugated );
        }
        else
        {
            // Handle 1x1
            Scale( d.GetLocal(jLoc,0), x1Loc );
        }
    }
}

template<typename F,typename FMain,Dist U,Dist V>
void
QuasiDiagonalScale
( LeftOrRight side, UpperOrLower uplo, 
  const AbstractDistMatrix<FMain>& d, const AbstractDistMatrix<F>& dSub, 
  DistMatrix<F,U,V>& X, bool conjugated )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiDiagonalScale"))
    const Grid& g = X.Grid();
    const Int colAlign = X.ColAlign();
    const Int rowAlign = X.RowAlign();
    if( side == LEFT )
    {
        const Int colStride = X.ColStride();
        DistMatrix<FMain,U,STAR> d_U_STAR(g);
        DistMatrix<F,U,STAR> dSub_U_STAR(g);
        d_U_STAR.AlignWith( X );
        dSub_U_STAR.AlignWith( X );
        d_U_STAR = d;
        dSub_U_STAR = dSub;
        if( colStride == 1 )
        {
            QuasiDiagonalScale
            ( side, uplo, 
              d_U_STAR.LockedMatrix(), dSub_U_STAR.LockedMatrix(),
              X.Matrix(), conjugated );
            return;
        }

        DistMatrix<FMain,U,STAR> dPrev_U_STAR(g), dNext_U_STAR(g);
        DistMatrix<F,U,STAR> dSubPrev_U_STAR(g), dSubNext_U_STAR(g);
        DistMatrix<F,U,V> XPrev(g), XNext(g);
        const Int colAlignPrev = Mod(colAlign+1,colStride);
        const Int colAlignNext = Mod(colAlign-1,colStride);
        dPrev_U_STAR.AlignCols( colAlignPrev );
        dNext_U_STAR.AlignCols( colAlignNext );
        dSubPrev_U_STAR.AlignCols( colAlignPrev );
        dSubNext_U_STAR.AlignCols( colAlignNext );
        XPrev.Align( colAlignPrev, rowAlign );
        XNext.Align( colAlignNext, rowAlign );
        dPrev_U_STAR = d;
        dNext_U_STAR = d;
        dSubPrev_U_STAR = dSub;
        dSubNext_U_STAR = dSub;
        XPrev = X;
        XNext = X;
        LeftQuasiDiagonalScale
        ( uplo, d_U_STAR, dPrev_U_STAR, dNext_U_STAR,
          dSub_U_STAR, dSubPrev_U_STAR, dSubNext_U_STAR,
          X, XPrev, XNext, conjugated );
    }
    else
    {
        const Int rowStride = X.RowStride();
        DistMatrix<FMain,V,STAR> d_V_STAR(g);
        DistMatrix<F,V,STAR> dSub_V_STAR(g);
        d_V_STAR.AlignWith( X );
        dSub_V_STAR.AlignWith( X );
        d_V_STAR = d;
        dSub_V_STAR = dSub;
        if( rowStride == 1 )
        {
            QuasiDiagonalScale
            ( side, uplo, 
              d_V_STAR.LockedMatrix(), dSub_V_STAR.LockedMatrix(),
              X.Matrix(), conjugated );
            return;
        }

        DistMatrix<FMain,V,STAR> dPrev_V_STAR(g), dNext_V_STAR(g);
        DistMatrix<F,V,STAR> dSubPrev_V_STAR(g), dSubNext_V_STAR(g);
        DistMatrix<F,U,V> XPrev(g), XNext(g);
        const Int rowAlignPrev = Mod(rowAlign+1,rowStride);
        const Int rowAlignNext = Mod(rowAlign-1,rowStride);
        dPrev_V_STAR.AlignCols( rowAlignPrev );
        dNext_V_STAR.AlignCols( rowAlignNext );
        dSubPrev_V_STAR.AlignCols( rowAlignPrev );
        dSubNext_V_STAR.AlignCols( rowAlignNext );
        XPrev.Align( colAlign, rowAlignPrev );
        XNext.Align( colAlign, rowAlignNext );
        dPrev_V_STAR = d;
        dNext_V_STAR = d;
        dSubPrev_V_STAR = dSub;
        dSubNext_V_STAR = dSub;
        XPrev = X;
        XNext = X;
        RightQuasiDiagonalScale
        ( uplo, d_V_STAR, dPrev_V_STAR, dNext_V_STAR,
          dSub_V_STAR, dSubPrev_V_STAR, dSubNext_V_STAR,
          X, XPrev, XNext, conjugated );
    }
}

#define PROTO_TYPES_DIST(F,FMain,U,V) \
  template void LeftQuasiDiagonalScale \
  ( UpperOrLower uplo, \
    const DistMatrix<FMain,U,STAR>& d, \
    const DistMatrix<FMain,U,STAR>& dPrev, \
    const DistMatrix<FMain,U,STAR>& dNext, \
    const DistMatrix<F,    U,STAR>& dSub, \
    const DistMatrix<F,    U,STAR>& dSubPrev, \
    const DistMatrix<F,    U,STAR>& dSubNext, \
          DistMatrix<F,U,V>& X, \
    const DistMatrix<F,U,V>& XPrev, \
    const DistMatrix<F,U,V>& XNext, \
          bool conjugated ); \
  template void RightQuasiDiagonalScale \
  ( UpperOrLower uplo, \
    const DistMatrix<FMain,V,STAR>& d, \
    const DistMatrix<FMain,V,STAR>& dPrev, \
    const DistMatrix<FMain,V,STAR>& dNext, \
    const DistMatrix<F,    V,STAR>& dSub, \
    const DistMatrix<F,    V,STAR>& dSubPrev, \
    const DistMatrix<F,    V,STAR>& dSubNext, \
          DistMatrix<F,U,V>& X, \
    const DistMatrix<F,U,V>& XPrev, \
    const DistMatrix<F,U,V>& XNext, \
          bool conjugated ); \
  template void QuasiDiagonalScale \
  ( LeftOrRight side, UpperOrLower uplo, \
    const AbstractDistMatrix<FMain>& d, const AbstractDistMatrix<F>& dSub, \
          DistMatrix<F,U,V>& X, bool conjugated );

#define PROTO_TYPES(F,FMain) \
  template void QuasiDiagonalScale \
  ( LeftOrRight side, UpperOrLower uplo, \
    const Matrix<FMain>& d, const Matrix<F>& dSub, \
          Matrix<F>& X, bool conjugated ); \
  PROTO_TYPES_DIST(F,FMain,MC,  MR  ) \
  PROTO_TYPES_DIST(F,FMain,MC,  STAR) \
  PROTO_TYPES_DIST(F,FMain,MD,  STAR) \
  PROTO_TYPES_DIST(F,FMain,MR,  MC  ) \
  PROTO_TYPES_DIST(F,FMain,MR,  STAR) \
  PROTO_TYPES_DIST(F,FMain,STAR,MC  ) \
  PROTO_TYPES_DIST(F,FMain,STAR,MD  ) \
  PROTO_TYPES_DIST(F,FMain,STAR,MR  ) \
  PROTO_TYPES_DIST(F,FMain,STAR,STAR) \
  PROTO_TYPES_DIST(F,FMain,STAR,VC  ) \
  PROTO_TYPES_DIST(F,FMain,STAR,VR  ) \
  PROTO_TYPES_DIST(F,FMain,VC,  STAR) \
  PROTO_TYPES_DIST(F,FMain,VR,  STAR) 

#define PROTO_REAL(F) PROTO_TYPES(F,F) 

#define PROTO_COMPLEX(F) \
  PROTO_TYPES(F,F) \
  PROTO_TYPES(F,Base<F>)

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
