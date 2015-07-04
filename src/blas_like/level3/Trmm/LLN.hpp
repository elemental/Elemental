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
#ifndef EL_TRMM_LLN_HPP
#define EL_TRMM_LLN_HPP

namespace El {
namespace trmm {

template<typename Ring>
inline void
LocalAccumulateLLN
( Orientation orientation, UnitOrNonUnit diag, Ring alpha,
  const DistMatrix<Ring,MC,  MR  >& L,
  const DistMatrix<Ring,STAR,MR  >& XTrans,
        DistMatrix<Ring,MC,  STAR>& Z )
{
    DEBUG_ONLY(
      CSE cse("trmm::LocalAccumulateLLN");
      AssertSameGrids( L, XTrans, Z );
      if( L.Height() != L.Width() ||
          L.Height() != XTrans.Width() ||
          L.Height() != Z.Height() ||
          XTrans.Height() != Z.Width() )
          LogicError
          ("Nonconformal: \n",DimsString(L,"L"),"\n",
           DimsString(XTrans,"X'"),"\n",DimsString(Z,"Z"));
      if( XTrans.RowAlign() != L.RowAlign() ||
          Z.ColAlign() != L.ColAlign() )
          LogicError("Partial matrix distributions are misaligned");
    )
    const Int m = Z.Height();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<Ring> D11(g);
    const Int ratio = Max( g.Height(), g.Width() );
    for( Int k=0; k<m; k+=ratio*bsize )
    {
        const Int nb = Min(ratio*bsize,m-k);

        auto L11 = L( IR(k,k+nb), IR(k,k+nb) );
        auto L21 = L( IR(k+nb,m), IR(k,k+nb) );

        auto X1Trans = XTrans( ALL, IR(k,k+nb) );

        auto Z1 = Z( IR(k,k+nb), ALL );
        auto Z2 = Z( IR(k+nb,m), ALL );

        D11.AlignWith( L11 );
        D11 = L11;
        MakeTrapezoidal( LOWER, D11 );
        if( diag == UNIT )
            FillDiagonal( D11, Ring(1) );
        LocalGemm( alpha, D11.N(), X1Trans.Orient(orientation), Ring(1), Z1 );
        LocalGemm( alpha, L21.N(), X1Trans.Orient(orientation), Ring(1), Z2 );
    }
}

template<typename Ring>
inline void
LLNA
( UnitOrNonUnit diag, 
  const AbstractDistMatrix<Ring>& LPre, AbstractDistMatrix<Ring>& XPre )
{
    DEBUG_ONLY(
      CSE cse("trmm::LLNA");
      AssertSameGrids( LPre, XPre );
      if( LPre.Height() != LPre.Width() || LPre.Width() != XPre.Height() )
          LogicError
          ("Nonconformal: \n",DimsString(LPre,"L"),"\n",DimsString(XPre,"X"))
    )
    const Int m = XPre.Height();
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    auto LPtr = ReadProxy<Ring,MC,MR>( &LPre );      auto& L = *LPtr;
    auto XPtr = ReadWriteProxy<Ring,MC,MR>( &XPre ); auto& X = *XPtr;

    DistMatrix<Ring,VR,  STAR> X1_VR_STAR(g);
    DistMatrix<Ring,STAR,MR  > X1Trans_STAR_MR(g);
    DistMatrix<Ring,MC,  STAR> Z1_MC_STAR(g);

    X1_VR_STAR.AlignWith( L );
    X1Trans_STAR_MR.AlignWith( L );
    Z1_MC_STAR.AlignWith( L );

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto X1 = X( ALL, IR(k,k+nb) );

        X1_VR_STAR = X1;
        Transpose( X1_VR_STAR, X1Trans_STAR_MR );
        Zeros( Z1_MC_STAR, m, nb );
        LocalAccumulateLLN
        ( TRANSPOSE, diag, Ring(1), L, X1Trans_STAR_MR, Z1_MC_STAR );
        Contract( Z1_MC_STAR, X1 );
    }
}

template<typename Ring>
inline void
LLNCOld
( UnitOrNonUnit diag, 
  const AbstractDistMatrix<Ring>& LPre, AbstractDistMatrix<Ring>& XPre )
{
    DEBUG_ONLY(
      CSE cse("trmm::LLNCOld");
      AssertSameGrids( LPre, XPre );
      if( LPre.Height() != LPre.Width() || LPre.Width() != XPre.Height() )
          LogicError
          ("Nonconformal: \n",DimsString(LPre,"L"),"\n",DimsString(XPre,"X"))
    )
    const Int m = XPre.Height();
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    auto LPtr = ReadProxy<Ring,MC,MR>( &LPre );      auto& L = *LPtr;
    auto XPtr = ReadWriteProxy<Ring,MC,MR>( &XPre ); auto& X = *XPtr;

    DistMatrix<Ring,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<Ring,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<Ring,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<Ring,MR,  STAR> D1Trans_MR_STAR(g);
    DistMatrix<Ring,MR,  MC  > D1Trans_MR_MC(g);
    DistMatrix<Ring,MC,  MR  > D1(g);

    const Int kLast = LastOffset( m, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto L10 = L( IR(k,k+nb), IR(0,k)    );
        auto L11 = L( IR(k,k+nb), IR(k,k+nb) );

        auto X0 = X( IR(0,k),    ALL );
        auto X1 = X( IR(k,k+nb), ALL );

        L11_STAR_STAR = L11;
        X1_STAR_VR = X1;
        LocalTrmm
        ( LEFT, LOWER, NORMAL, diag, Ring(1), L11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;

        L10_STAR_MC.AlignWith( X0 );
        L10_STAR_MC = L10;
        D1Trans_MR_STAR.AlignWith( X1 );
        LocalGemm( Ring(1), X0.T(), L10_STAR_MC.T(), D1Trans_MR_STAR );
        D1Trans_MR_MC.AlignWith( X1 );
        Contract( D1Trans_MR_STAR, D1Trans_MR_MC );
        D1.AlignWith( X1 );
        Zeros( D1, nb, n );
        Transpose( D1Trans_MR_MC.Matrix(), D1.Matrix() );
        X1 += D1;
    }
}

template<typename Ring>
inline void
LLNC
( UnitOrNonUnit diag, 
  const AbstractDistMatrix<Ring>& LPre, AbstractDistMatrix<Ring>& XPre )
{
    DEBUG_ONLY(
      CSE cse("trmm::LLNC");
      AssertSameGrids( LPre, XPre );
      if( LPre.Height() != LPre.Width() || LPre.Width() != XPre.Height() )
          LogicError
          ("Nonconformal: \n",DimsString(LPre,"L"),"\n",DimsString(XPre,"X"))
    )
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    auto LPtr = ReadProxy<Ring,MC,MR>( &LPre );      auto& L = *LPtr;
    auto XPtr = ReadWriteProxy<Ring,MC,MR>( &XPre ); auto& X = *XPtr;

    DistMatrix<Ring,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<Ring,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<Ring,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<Ring,MR,  STAR> X1Trans_MR_STAR(g);

    const Int kLast = LastOffset( m, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto L11 = L( IR(k,k+nb), IR(k,k+nb) );
        auto L21 = L( IR(k+nb,m), IR(k,k+nb) );

        auto X1 = X( IR(k,k+nb), ALL );
        auto X2 = X( IR(k+nb,m), ALL );

        L21_MC_STAR.AlignWith( X2 );
        L21_MC_STAR = L21;
        X1Trans_MR_STAR.AlignWith( X2 );
        Transpose( X1, X1Trans_MR_STAR );
        LocalGemm( Ring(1), L21_MC_STAR.N(), X1Trans_MR_STAR.T(), Ring(1), X2 );

        L11_STAR_STAR = L11;
        X1_STAR_VR.AlignWith( X1 );
        Transpose( X1Trans_MR_STAR, X1_STAR_VR );
        LocalTrmm
        ( LEFT, LOWER, NORMAL, diag, Ring(1), L11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
    }
}

// Left Lower Normal (Non)Unit Trmm 
//   X := tril(L)  X, or
//   X := trilu(L) X
template<typename Ring>
inline void
LLN
( UnitOrNonUnit diag, 
  const AbstractDistMatrix<Ring>& L, AbstractDistMatrix<Ring>& X )
{
    DEBUG_ONLY(CSE cse("trmm::LLN"))
    // TODO: Come up with a better routing mechanism
    if( L.Height() > 5*X.Width() )
        LLNA( diag, L, X );
    else
        LLNC( diag, L, X );
}

} // namespace trmm
} // namespace El

#endif // ifndef EL_TRMM_LLN_HPP
