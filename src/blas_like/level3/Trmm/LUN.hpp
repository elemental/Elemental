/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2013, The University of Texas at Austin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRMM_LUN_HPP
#define EL_TRMM_LUN_HPP

namespace El {
namespace trmm {

template<typename T>
inline void
LocalAccumulateLUN
( Orientation orientation,
  UnitOrNonUnit diag,
  T alpha,
  const DistMatrix<T,MC,  MR  >& U,
  const DistMatrix<T,STAR,MR  >& XTrans,
        DistMatrix<T,MC,  STAR>& Z )
{
    DEBUG_ONLY(
      CSE cse("trmm::LocalAccumulateLUN");
      AssertSameGrids( U, XTrans, Z );
      if( U.Height() != U.Width() ||
          U.Height() != XTrans.Width() ||
          U.Height() != Z.Height() ||
          XTrans.Height() != Z.Width() )
          LogicError
          ("Nonconformal:\n",
           DimsString(U,"U"),"\n",
           DimsString(XTrans,"X'[* ,MR]"),"\n",
           DimsString(Z,"Z[MC,* ]"));
      if( XTrans.RowAlign() != U.RowAlign() ||
          Z.ColAlign() != U.ColAlign() )
          LogicError("Partial matrix distributions are misaligned");
    )
    const Int m = Z.Height();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<T> D11(g);

    const Int ratio = Max( g.Height(), g.Width() );
    for( Int k=0; k<m; k+=ratio*bsize )
    {
        const Int nb = Min(ratio*bsize,m-k);

        auto U01 = U( IR(0,k),    IR(k,k+nb) );
        auto U11 = U( IR(k,k+nb), IR(k,k+nb) );

        auto X1Trans = XTrans( ALL, IR(k,k+nb) );

        auto Z0 = Z( IR(0,k),    ALL );
        auto Z1 = Z( IR(k,k+nb), ALL );

        D11.AlignWith( U11 );
        D11 = U11;
        MakeTrapezoidal( UPPER, D11 );
        if( diag == UNIT )
            FillDiagonal( D11, T(1) );
        LocalGemm( NORMAL, orientation, alpha, D11, X1Trans, T(1), Z1 );
        LocalGemm( NORMAL, orientation, alpha, U01, X1Trans, T(1), Z0 );
    }
}

template<typename T>
inline void
LUNA
( UnitOrNonUnit diag, 
  const ElementalMatrix<T>& UPre,
        ElementalMatrix<T>& XPre )
{
    DEBUG_ONLY(
      CSE cse("trmm::LUNA");
      AssertSameGrids( UPre, XPre );
      if( UPre.Height() != UPre.Width() || UPre.Width() != XPre.Height() )
          LogicError
          ("Nonconformal:\n",DimsString(UPre,"U"),"\n",DimsString(XPre,"X"));
    )
    const Int m = XPre.Height();
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();

    DistMatrixReadProxy<T,T,MC,MR> UProx( UPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<T,VR,  STAR> X1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > X1Trans_STAR_MR(g);
    DistMatrix<T,MC,  STAR> Z1_MC_STAR(g);

    X1_VR_STAR.AlignWith( U );
    X1Trans_STAR_MR.AlignWith( U );
    Z1_MC_STAR.AlignWith( U );

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto X1 = X( ALL, IR(k,k+nb) );

        X1_VR_STAR = X1;
        Transpose( X1_VR_STAR, X1Trans_STAR_MR );
        Zeros( Z1_MC_STAR, m, nb );
        LocalAccumulateLUN
        ( TRANSPOSE, diag, T(1), U, X1Trans_STAR_MR, Z1_MC_STAR );

        Contract( Z1_MC_STAR, X1 );
    }
}

template<typename T>
inline void
LUNCOld
( UnitOrNonUnit diag, 
  const ElementalMatrix<T>& UPre,
        ElementalMatrix<T>& XPre )
{
    DEBUG_ONLY(
      CSE cse("trmm::LUNCOld");
      AssertSameGrids( UPre, XPre );
      if( UPre.Height() != UPre.Width() || UPre.Width() != XPre.Height() )
          LogicError
          ("Nonconformal:\n",DimsString(UPre,"U"),"\n",DimsString(XPre,"X"));
    )
    const Int m = XPre.Height();
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();

    DistMatrixReadProxy<T,T,MC,MR> UProx( UPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<T,STAR,MC  > U12_STAR_MC(g);
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> D1Trans_MR_STAR(g);
    DistMatrix<T,MR,  MC  > D1Trans_MR_MC(g);
    DistMatrix<T,MC,  MR  > D1(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto U11 = U( IR(k,k+nb), IR(k,k+nb) );
        auto U12 = U( IR(k,k+nb), IR(k+nb,m) );

        auto X1 = X( IR(k,k+nb), ALL );
        auto X2 = X( IR(k+nb,m), ALL );

        X1_STAR_VR = X1;
        U11_STAR_STAR = U11;
        LocalTrmm
        ( LEFT, UPPER, NORMAL, diag, T(1), U11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
 
        U12_STAR_MC.AlignWith( X2 );
        U12_STAR_MC = U12;
        D1Trans_MR_STAR.AlignWith( X1 );
        LocalGemm
        ( TRANSPOSE, TRANSPOSE, T(1), X2, U12_STAR_MC, D1Trans_MR_STAR );
        D1Trans_MR_MC.AlignWith( X1 );
        Contract( D1Trans_MR_STAR, D1Trans_MR_MC );
        D1.AlignWith( X1 );
        Zeros( D1, nb, n );
        Transpose( D1Trans_MR_MC.Matrix(), D1.Matrix() );
        Axpy( T(1), D1, X1 );
    }
}

template<typename T>
inline void
LUNC
( UnitOrNonUnit diag, 
  const ElementalMatrix<T>& UPre,
        ElementalMatrix<T>& XPre )
{
    DEBUG_ONLY(
      CSE cse("trmm::LUNC");
      AssertSameGrids( UPre, XPre );
      if( UPre.Height() != UPre.Width() || UPre.Width() != XPre.Height() )
          LogicError
          ("Nonconformal:\n",DimsString(UPre,"U"),"\n",DimsString(XPre,"X"));
    )
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();

    DistMatrixReadProxy<T,T,MC,MR> UProx( UPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<T,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> X1Trans_MR_STAR(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto U01 = U( IR(0,k),    IR(k,k+nb) );
        auto U11 = U( IR(k,k+nb), IR(k,k+nb) );

        auto X0 = X( IR(0,k),    ALL );
        auto X1 = X( IR(k,k+nb), ALL );

        U01_MC_STAR.AlignWith( X0 );
        U01_MC_STAR = U01;
        X1Trans_MR_STAR.AlignWith( X0 );
        Transpose( X1, X1Trans_MR_STAR );
        LocalGemm
        ( NORMAL, TRANSPOSE, T(1), U01_MC_STAR, X1Trans_MR_STAR, T(1), X0 );

        U11_STAR_STAR = U11;
        X1_STAR_VR.AlignWith( X1 );
        Transpose( X1Trans_MR_STAR, X1_STAR_VR );
        LocalTrmm( LEFT, UPPER, NORMAL, diag, T(1), U11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
    }
}

// Left Upper Normal (Non)Unit Trmm
//   X := triu(U)  X, or
//   X := triuu(U) X
template<typename T>
inline void
LUN
( UnitOrNonUnit diag,
  const ElementalMatrix<T>& U,
        ElementalMatrix<T>& X )
{
    DEBUG_ONLY(CSE cse("trmm::LUN"))
    // TODO: Come up with a better routing mechanism
    if( U.Height() > 5*X.Width() )
        LUNA( diag, U, X );
    else
        LUNC( diag, U, X );
}

} // namespace trmm
} // namespace El

#endif // ifndef EL_TRMM_LUN_HPP
