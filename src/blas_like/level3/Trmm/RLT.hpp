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
#ifndef EL_TRMM_RLT_HPP
#define EL_TRMM_RLT_HPP

namespace El {
namespace trmm {

template<typename Ring>
inline void
LocalAccumulateRLT
( UnitOrNonUnit diag, Ring alpha,
  const DistMatrix<Ring>& L,
  const DistMatrix<Ring,MR,STAR>& XTrans,
        DistMatrix<Ring,MC,STAR>& ZTrans )
{
    DEBUG_ONLY(
      CSE cse("trmm::LocalAccumulateRLT");
      AssertSameGrids( L, XTrans, ZTrans );
      if( L.Height() != L.Width() ||
          L.Height() != XTrans.Height() ||
          L.Height() != ZTrans.Height() ||
          XTrans.Width() != ZTrans.Width() )
          LogicError
          ("Nonconformal:\n",DimsString(L,"L"),"\n",
           DimsString(XTrans,"X'"),"\n",DimsString(ZTrans,"Z'"));
      if( XTrans.ColAlign() != L.RowAlign() ||
          ZTrans.ColAlign() != L.ColAlign() )
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

        auto X1Trans = XTrans( IR(k,k+nb), ALL );

        auto Z1Trans = ZTrans( IR(k,k+nb), ALL );
        auto Z2Trans = ZTrans( IR(k+nb,m), ALL );

        D11.AlignWith( L11 );
        D11 = L11;
        MakeTrapezoidal( LOWER, D11 );
        if( diag == UNIT )
            FillDiagonal( D11, Ring(1) );
        LocalGemm( alpha, D11.N(), X1Trans.N(), Ring(1), Z1Trans );
        LocalGemm( alpha, L21.N(), X1Trans.N(), Ring(1), Z2Trans );
    }
}

template<typename Ring>
inline void
RLTA
( Orientation orientation, UnitOrNonUnit diag,
  const AbstractDistMatrix<Ring>& LPre, AbstractDistMatrix<Ring>& XPre )
{
    DEBUG_ONLY(
      CSE cse("trmm::RLTA");
      AssertSameGrids( LPre, XPre );
      // TODO: More error checks
    )
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();
    const bool conjugate = ( orientation == ADJOINT );

    auto LPtr = ReadProxy<Ring,MC,MR>( &LPre );      auto& L = *LPtr;
    auto XPtr = ReadWriteProxy<Ring,MC,MR>( &XPre ); auto& X = *XPtr;

    DistMatrix<Ring,MR,  STAR> X1Trans_MR_STAR(g);
    DistMatrix<Ring,MC,  STAR> Z1Trans_MC_STAR(g);
    DistMatrix<Ring,MC,  MR  > Z1Trans(g);
    DistMatrix<Ring,MR,  MC  > Z1Trans_MR_MC(g);

    X1Trans_MR_STAR.AlignWith( L );
    Z1Trans_MC_STAR.AlignWith( L );

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto X1 = X( IR(k,k+nb), ALL );

        Transpose( X1, X1Trans_MR_STAR, conjugate );
        Zeros( Z1Trans_MC_STAR, X1.Width(), X1.Height() );
        LocalAccumulateRLT
        ( diag, Ring(1), L, X1Trans_MR_STAR, Z1Trans_MC_STAR );

        Contract( Z1Trans_MC_STAR, Z1Trans );
        Z1Trans_MR_MC.AlignWith( X1 );
        Z1Trans_MR_MC = Z1Trans;
        Transpose( Z1Trans_MR_MC.Matrix(), X1.Matrix(), conjugate );
    }
}

template<typename Ring>
inline void
RLTC
( Orientation orientation, UnitOrNonUnit diag,
  const AbstractDistMatrix<Ring>& LPre, AbstractDistMatrix<Ring>& XPre )
{
    DEBUG_ONLY(
      CSE cse("trmm::RLTC");
      AssertSameGrids( LPre, XPre );
      if( orientation == NORMAL )
          LogicError("Expected Adjoint/Transpose option");
      if( LPre.Height() != LPre.Width() || XPre.Width() != LPre.Height() )
          LogicError
          ("Nonconformal: \n",DimsString(LPre,"L"),"\n",DimsString(XPre,"X"))
    )
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();
    const bool conjugate = ( orientation == ADJOINT );

    auto LPtr = ReadProxy<Ring,MC,MR>( &LPre );      auto& L = *LPtr;
    auto XPtr = ReadWriteProxy<Ring,MC,MR>( &XPre ); auto& X = *XPtr;

    DistMatrix<Ring,MR,  STAR> L10Trans_MR_STAR(g);
    DistMatrix<Ring,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<Ring,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<Ring,MC,  STAR> D1_MC_STAR(g);

    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);
         
        auto L10 = L( IR(k,k+nb), IR(0,k)    );
        auto L11 = L( IR(k,k+nb), IR(k,k+nb) );

        auto X0 = X( ALL, IR(0,k)    );
        auto X1 = X( ALL, IR(k,k+nb) );

        X1_VC_STAR = X1;
        L11_STAR_STAR = L11;
        LocalTrmm
        ( RIGHT, LOWER, orientation, diag, Ring(1), L11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
 
        L10Trans_MR_STAR.AlignWith( X0 );
        Transpose( L10, L10Trans_MR_STAR, conjugate );
        D1_MC_STAR.AlignWith( X1 );
        LocalGemm( Ring(1), X0.N(), L10Trans_MR_STAR.N(), D1_MC_STAR );
        AxpyContract( Ring(1), D1_MC_STAR, X1 );
    }
}

// Right Lower Adjoint/Transpose (Non)Unit Trmm
//   X := X tril(L)^T,
//   X := X tril(L)^H,
//   X := X trilu(L)^T, or
//   X := X trilu(L)^H
template<typename Ring>
inline void
RLT
( Orientation orientation, UnitOrNonUnit diag,
  const AbstractDistMatrix<Ring>& L, AbstractDistMatrix<Ring>& X )
{
    DEBUG_ONLY(CSE cse("trmm::RLT"))
    // TODO: Come up with a better routing mechanism
    if( L.Height() > 5*X.Height() )
        RLTA( orientation, diag, L, X );
    else
        RLTC( orientation, diag, L, X );
}

} // namespace trmm
} // namespace El

#endif // ifndef EL_TRMM_RLT_HPP
