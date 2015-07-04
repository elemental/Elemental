/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   Copyright (c) 2013, The University of Texas at Austin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_TRMM_RLN_HPP
#define EL_TRMM_RLN_HPP

namespace El {
namespace trmm {

template<typename Ring>
inline void
LocalAccumulateRLN
( Orientation orient, UnitOrNonUnit diag, Ring alpha,
  const DistMatrix<Ring,MC,  MR  >& L,
  const DistMatrix<Ring,STAR,MC  >& X,
        DistMatrix<Ring,MR,  STAR>& ZTrans )
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

    DistMatrix<Ring> D11(g);

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
            FillDiagonal( D11, Ring(1) );
        LocalGemm
        ( alpha, D11.Orient(orient), X1.Orient(orient), Ring(1), Z1Trans );
        LocalGemm
        ( alpha, L21.Orient(orient), X2.Orient(orient), Ring(1), Z1Trans );
    }
}

template<typename Ring>
inline void
RLNA
( UnitOrNonUnit diag, 
  const AbstractDistMatrix<Ring>& LPre, AbstractDistMatrix<Ring>& XPre )
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

    auto LPtr = ReadProxy<Ring,MC,MR>( &LPre );      auto& L = *LPtr;
    auto XPtr = ReadWriteProxy<Ring,MC,MR>( &XPre ); auto& X = *XPtr;

    DistMatrix<Ring,STAR,VC  > X1_STAR_VC(g);
    DistMatrix<Ring,STAR,MC  > X1_STAR_MC(g);
    DistMatrix<Ring,MR,  STAR> Z1Trans_MR_STAR(g);
    DistMatrix<Ring,MR,  MC  > Z1Trans_MR_MC(g);

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
        ( TRANSPOSE, diag, Ring(1), L, X1_STAR_MC, Z1Trans_MR_STAR );

        Z1Trans_MR_MC.AlignWith( X1 );
        Contract( Z1Trans_MR_STAR, Z1Trans_MR_MC );
        Transpose( Z1Trans_MR_MC.Matrix(), X1.Matrix() );
    }
}

template<typename Ring>
inline void
RLNCOld
( UnitOrNonUnit diag, 
  const AbstractDistMatrix<Ring>& LPre, AbstractDistMatrix<Ring>& XPre )
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
 
    auto LPtr = ReadProxy<Ring,MC,MR>( &LPre );      auto& L = *LPtr;
    auto XPtr = ReadWriteProxy<Ring,MC,MR>( &XPre ); auto& X = *XPtr;

    DistMatrix<Ring,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<Ring,MR,  STAR> L21_MR_STAR(g);
    DistMatrix<Ring,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<Ring,MC,  STAR> D1_MC_STAR(g);

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
        ( RIGHT, LOWER, NORMAL, diag, Ring(1), L11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
 
        L21_MR_STAR.AlignWith( X2 );
        L21_MR_STAR = L21;
        D1_MC_STAR.AlignWith( X1 );
        LocalGemm( Ring(1), X2.N(), L21_MR_STAR.N(), D1_MC_STAR );
        AxpyContract( Ring(1), D1_MC_STAR, X1 );
    }
}

template<typename Ring>
inline void
RLNC
( UnitOrNonUnit diag, 
  const AbstractDistMatrix<Ring>& LPre, AbstractDistMatrix<Ring>& XPre )
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

    auto LPtr = ReadProxy<Ring,MC,MR>( &LPre );      auto& L = *LPtr;
    auto XPtr = ReadWriteProxy<Ring,MC,MR>( &XPre ); auto& X = *XPtr;

    DistMatrix<Ring,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<Ring,MR,  STAR> L10Trans_MR_STAR(g);
    DistMatrix<Ring,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<Ring,MC,  STAR> X1_MC_STAR(g);

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
        ( Ring(1), X1_MC_STAR.N(), L10Trans_MR_STAR.T(), Ring(1), X0 );

        L11_STAR_STAR = L11;
        X1_VC_STAR.AlignWith( X1 );
        X1_VC_STAR = X1_MC_STAR;
        LocalTrmm
        ( RIGHT, LOWER, NORMAL, diag, Ring(1), L11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
    }
}

// Right Lower Normal (Non)Unit Trmm
//   X := X tril(L), and
//   X := X trilu(L)
template<typename Ring>
inline void
RLN
( UnitOrNonUnit diag, 
  const AbstractDistMatrix<Ring>& L, AbstractDistMatrix<Ring>& X )
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
