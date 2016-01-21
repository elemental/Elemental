/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2013, The University of Texas at Austin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRMM_RLN_HPP
#define EL_TRMM_RLN_HPP

namespace El {
namespace trmm {

template<typename T>
inline void
LocalAccumulateRLN
( Orientation orientation,
  UnitOrNonUnit diag,
  T alpha,
  const DistMatrix<T,MC,  MR  >& L,
  const DistMatrix<T,STAR,MC  >& X,
        DistMatrix<T,MR,  STAR>& ZTrans )
{
    DEBUG_ONLY(
      CSE cse("trmm::LocalAccumulateRLN");
      AssertSameGrids( L, X, ZTrans );
      if( L.Height() != L.Width() ||
          L.Height() != X.Width() ||
          L.Height() != ZTrans.Height() )
          LogicError
          ("Nonconformal:\n",
           DimsString(L,"L"),"\n",
           DimsString(X,"X[* ,MC]"),"\n",
           DimsString(ZTrans,"Z'[MR,* ]"));
      if( X.RowAlign() != L.ColAlign() ||
          ZTrans.ColAlign() != L.RowAlign() )
          LogicError("Partial matrix distributions are misaligned");
    )
    const Int m = ZTrans.Height();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<T> D11(g);

    const Int ratio = Max( g.Height(), g.Width() );
    for( Int k=0; k<m; k+=ratio*bsize )
    {
        const Int nb = Min(ratio*bsize,m-k);

        auto L11 = L( IR(k,k+nb), IR(k,k+nb) );
        auto L21 = L( IR(k+nb,m), IR(k,k+nb) );

        auto X1 = X( ALL, IR(k,k+nb) );
        auto X2 = X( ALL, IR(k+nb,m) );

        auto Z1Trans = ZTrans( IR(k,k+nb), ALL );

        D11.AlignWith( L11 );
        D11 = L11;
        MakeTrapezoidal( LOWER, D11 );
        if( diag == UNIT )
            FillDiagonal( D11, T(1) );
        LocalGemm( orientation, orientation, alpha, D11, X1, T(1), Z1Trans );
        LocalGemm( orientation, orientation, alpha, L21, X2, T(1), Z1Trans );
    }
}

template<typename T>
inline void
RLNA
( UnitOrNonUnit diag, 
  const ElementalMatrix<T>& LPre,
        ElementalMatrix<T>& XPre )
{
    DEBUG_ONLY(
      CSE cse("trmm::RLNA");
      AssertSameGrids( LPre, XPre );
      // TODO: More checks
    )
    const Int m = XPre.Height();
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    DistMatrixReadProxy<T,T,MC,MR> LProx( LPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> XProx( XPre );
    auto& L = LProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<T,STAR,VC  > X1_STAR_VC(g);
    DistMatrix<T,STAR,MC  > X1_STAR_MC(g);
    DistMatrix<T,MR,  STAR> Z1Trans_MR_STAR(g);
    DistMatrix<T,MR,  MC  > Z1Trans_MR_MC(g);

    X1_STAR_VC.AlignWith( L );
    X1_STAR_MC.AlignWith( L );
    Z1Trans_MR_STAR.AlignWith( L );

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto X1 = X( IR(k,k+nb), ALL );

        X1_STAR_VC = X1;
        X1_STAR_MC = X1_STAR_VC;

        Zeros( Z1Trans_MR_STAR, n, nb );
        LocalAccumulateRLN
        ( TRANSPOSE, diag, T(1), L, X1_STAR_MC, Z1Trans_MR_STAR );

        Z1Trans_MR_MC.AlignWith( X1 );
        Contract( Z1Trans_MR_STAR, Z1Trans_MR_MC );
        Transpose( Z1Trans_MR_MC.Matrix(), X1.Matrix() );
    }
}

template<typename T>
inline void
RLNCOld
( UnitOrNonUnit diag, 
  const ElementalMatrix<T>& LPre,
        ElementalMatrix<T>& XPre )
{
    DEBUG_ONLY(
      CSE cse("trmm::RLNCOld");
      AssertSameGrids( LPre, XPre );
      if( LPre.Height() != LPre.Width() || XPre.Width() != LPre.Height() )
          LogicError
          ("Nonconformal:\n",DimsString(LPre,"L"),"\n",DimsString(XPre,"X"));
    )
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    DistMatrixReadProxy<T,T,MC,MR> LProx( LPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> XProx( XPre );
    auto& L = LProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,MR,  STAR> L21_MR_STAR(g);
    DistMatrix<T,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<T,MC,  STAR> D1_MC_STAR(g);

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto L11 = L( IR(k,k+nb), IR(k,k+nb) );
        auto L21 = L( IR(k+nb,n), IR(k,k+nb) );

        auto X1 = X( ALL, IR(k,k+nb) );
        auto X2 = X( ALL, IR(k+nb,n) );

        X1_VC_STAR = X1;
        L11_STAR_STAR = L11;
        LocalTrmm
        ( RIGHT, LOWER, NORMAL, diag, T(1), L11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
 
        L21_MR_STAR.AlignWith( X2 );
        L21_MR_STAR = L21;
        D1_MC_STAR.AlignWith( X1 );
        LocalGemm( NORMAL, NORMAL, T(1), X2, L21_MR_STAR, D1_MC_STAR );
        AxpyContract( T(1), D1_MC_STAR, X1 );
    }
}

template<typename T>
inline void
RLNC
( UnitOrNonUnit diag, 
  const ElementalMatrix<T>& LPre,
        ElementalMatrix<T>& XPre )
{
    DEBUG_ONLY(
      CSE cse("trmm::RLNC");
      AssertSameGrids( LPre, XPre );
      if( LPre.Height() != LPre.Width() || XPre.Width() != LPre.Height() )
          LogicError
          ("Nonconformal:\n",DimsString(LPre,"L"),"\n",DimsString(XPre,"X"));
    )
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    DistMatrixReadProxy<T,T,MC,MR> LProx( LPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> XProx( XPre );
    auto& L = LProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,MR,  STAR> L10Trans_MR_STAR(g);
    DistMatrix<T,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<T,MC,  STAR> X1_MC_STAR(g);

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto L10 = L( IR(k,k+nb), IR(0,k)    );
        auto L11 = L( IR(k,k+nb), IR(k,k+nb) );

        auto X0 = X( ALL, IR(0,k)    );
        auto X1 = X( ALL, IR(k,k+nb) );

        X1_MC_STAR.AlignWith( X0 );
        X1_MC_STAR = X1;
        L10Trans_MR_STAR.AlignWith( X0 );
        Transpose( L10, L10Trans_MR_STAR );
        LocalGemm
        ( NORMAL, TRANSPOSE, T(1), X1_MC_STAR, L10Trans_MR_STAR, T(1), X0 );

        L11_STAR_STAR = L11;
        X1_VC_STAR.AlignWith( X1 );
        X1_VC_STAR = X1_MC_STAR;
        LocalTrmm
        ( RIGHT, LOWER, NORMAL, diag, T(1), L11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
    }
}

// Right Lower Normal (Non)Unit Trmm
//   X := X tril(L), and
//   X := X trilu(L)
template<typename T>
inline void
RLN
( UnitOrNonUnit diag, 
  const ElementalMatrix<T>& L,
        ElementalMatrix<T>& X )
{
    DEBUG_ONLY(CSE cse("trmm::RLN"))
    // TODO: Come up with a better routing mechanism
    if( L.Height() > 5*X.Height() )
        RLNA( diag, L, X );
    else
        RLNC( diag, L, X );
}

} // namespace trmm
} // namespace El

#endif // ifndef EL_TRMM_RLN_HPP
