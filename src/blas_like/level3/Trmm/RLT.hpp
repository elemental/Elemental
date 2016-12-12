/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2013, The University of Texas at Austin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRMM_RLT_HPP
#define EL_TRMM_RLT_HPP

namespace El {
namespace trmm {

template<typename T>
void LocalAccumulateRLT
( UnitOrNonUnit diag,
  T alpha,
  const DistMatrix<T>& L,
  const DistMatrix<T,MR,STAR>& XTrans,
        DistMatrix<T,MC,STAR>& ZTrans )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
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

    DistMatrix<T> D11(g);

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
            FillDiagonal( D11, T(1) );
        LocalGemm( NORMAL, NORMAL, alpha, D11, X1Trans, T(1), Z1Trans );
        LocalGemm( NORMAL, NORMAL, alpha, L21, X1Trans, T(1), Z2Trans );
    }
}

template<typename T>
void RLTA
( Orientation orientation,
  UnitOrNonUnit diag,
  const AbstractDistMatrix<T>& LPre,
        AbstractDistMatrix<T>& XPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( LPre, XPre );
      // TODO: More error checks
    )
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();
    const bool conjugate = ( orientation == ADJOINT );

    DistMatrixReadProxy<T,T,MC,MR> LProx( LPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> XProx( XPre );
    auto& L = LProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<T,MR,  STAR> X1Trans_MR_STAR(g);
    DistMatrix<T,MC,  STAR> Z1Trans_MC_STAR(g);
    DistMatrix<T,MC,  MR  > Z1Trans(g);
    DistMatrix<T,MR,  MC  > Z1Trans_MR_MC(g);

    X1Trans_MR_STAR.AlignWith( L );
    Z1Trans_MC_STAR.AlignWith( L );

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto X1 = X( IR(k,k+nb), ALL );

        Transpose( X1, X1Trans_MR_STAR, conjugate );
        Z1Trans_MC_STAR.Resize( X1.Width(), X1.Height() );
        Zero( Z1Trans_MC_STAR );
        LocalAccumulateRLT
        ( diag, T(1), L, X1Trans_MR_STAR, Z1Trans_MC_STAR );

        Contract( Z1Trans_MC_STAR, Z1Trans );
        Z1Trans_MR_MC.AlignWith( X1 );
        Z1Trans_MR_MC = Z1Trans;
        Transpose( Z1Trans_MR_MC.Matrix(), X1.Matrix(), conjugate );
    }
}

template<typename T>
void RLTC
( Orientation orientation,
  UnitOrNonUnit diag,
  const AbstractDistMatrix<T>& LPre,
        AbstractDistMatrix<T>& XPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
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

    DistMatrixReadProxy<T,T,MC,MR> LProx( LPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> XProx( XPre );
    auto& L = LProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<T,MR,  STAR> L10Trans_MR_STAR(g);
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<T,MC,  STAR> D1_MC_STAR(g);

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
        ( RIGHT, LOWER, orientation, diag, T(1), L11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
 
        L10Trans_MR_STAR.AlignWith( X0 );
        Transpose( L10, L10Trans_MR_STAR, conjugate );
        D1_MC_STAR.AlignWith( X1 );
        LocalGemm( NORMAL, NORMAL, T(1), X0, L10Trans_MR_STAR, D1_MC_STAR );
        AxpyContract( T(1), D1_MC_STAR, X1 );
    }
}

// Right Lower Adjoint/Transpose (Non)Unit Trmm
//   X := X tril(L)^T,
//   X := X tril(L)^H,
//   X := X trilu(L)^T, or
//   X := X trilu(L)^H
template<typename T>
void RLT
( Orientation orientation,
  UnitOrNonUnit diag,
  const AbstractDistMatrix<T>& L,
        AbstractDistMatrix<T>& X )
{
    EL_DEBUG_CSE
    // TODO: Come up with a better routing mechanism
    if( L.Height() > 5*X.Height() )
        RLTA( orientation, diag, L, X );
    else
        RLTC( orientation, diag, L, X );
}

} // namespace trmm
} // namespace El

#endif // ifndef EL_TRMM_RLT_HPP
