/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2013, The University of Texas at Austin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRMM_RUT_HPP
#define EL_TRMM_RUT_HPP

namespace El {
namespace trmm {

template<typename T>
void LocalAccumulateRUT
( UnitOrNonUnit diag,
  T alpha,
  const DistMatrix<T>& U,
  const DistMatrix<T,MR,STAR>& XTrans,
        DistMatrix<T,MC,STAR>& ZTrans )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( U, XTrans, ZTrans );
      if( U.Height() != U.Width() ||
          U.Height() != XTrans.Height() ||
          U.Height() != ZTrans.Height() ||
          XTrans.Width() != ZTrans.Width() )
          LogicError
          ("Nonconformal: \n",DimsString(U,"U"),"\n",
           DimsString(XTrans,"X'"),"\n",DimsString(ZTrans,"Z'"));
      if( XTrans.ColAlign() != U.RowAlign() ||
          ZTrans.ColAlign() != U.ColAlign() )
          LogicError("Partial matrix distributions are misaligned");
    )
    const Int m = ZTrans.Height();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<T> D11(g);

    const Int ratio = Max( g.Height(), g.Width() );
    for( Int k=0; k<m; k+=ratio*bsize )
    {
        const Int nb = Min(ratio*bsize,m-k);

        auto U01 = U( IR(0,k),    IR(k,k+nb) );
        auto U11 = U( IR(k,k+nb), IR(k,k+nb) );

        auto X1Trans = XTrans( IR(k,k+nb), ALL );
   
        auto Z0Trans = ZTrans( IR(0,k),    ALL );
        auto Z1Trans = ZTrans( IR(k,k+nb), ALL );

        D11.AlignWith( U11 );
        D11 = U11;
        MakeTrapezoidal( UPPER, D11 );
        if( diag == UNIT )
            FillDiagonal( D11, T(1) );
        LocalGemm( NORMAL, NORMAL, alpha, D11, X1Trans, T(1), Z1Trans );
        LocalGemm( NORMAL, NORMAL, alpha, U01, X1Trans, T(1), Z0Trans );
    }
}

template<typename T>
void RUTA
( Orientation orientation,
  UnitOrNonUnit diag,
  const AbstractDistMatrix<T>& UPre,
        AbstractDistMatrix<T>& XPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( UPre, XPre );
      // TODO: More input checks
    )
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();
    const bool conjugate = ( orientation == ADJOINT );

    DistMatrixReadProxy<T,T,MC,MR> UProx( UPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<T,MR,  STAR> X1Trans_MR_STAR(g);
    DistMatrix<T,MC,  STAR> Z1Trans_MC_STAR(g);
    DistMatrix<T,MC,  MR  > Z1Trans(g);
    DistMatrix<T,MR,  MC  > Z1Trans_MR_MC(g);

    X1Trans_MR_STAR.AlignWith( U );
    Z1Trans_MC_STAR.AlignWith( U );

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto X1 = X( IR(k,k+nb), ALL );

        Transpose( X1, X1Trans_MR_STAR, conjugate );
        Z1Trans_MC_STAR.Resize( X.Width(), nb );
        Zero( Z1Trans_MC_STAR );
        LocalAccumulateRUT
        ( diag, T(1), U, X1Trans_MR_STAR, Z1Trans_MC_STAR );

        Contract( Z1Trans_MC_STAR, Z1Trans );
        Z1Trans_MR_MC.AlignWith( X1 );
        Z1Trans_MR_MC = Z1Trans;
        Transpose( Z1Trans_MR_MC.Matrix(), X1.Matrix(), conjugate );
    }
}

template<typename T>
void RUTC
( Orientation orientation,
  UnitOrNonUnit diag,
  const AbstractDistMatrix<T>& UPre,
        AbstractDistMatrix<T>& XPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( UPre, XPre );
      if( orientation == NORMAL )
          LogicError("Expected Adjoint/Transpose option");
      if( UPre.Height() != UPre.Width() || XPre.Width() != UPre.Height() )
          LogicError
          ("Nonconformal: \n",DimsString(UPre,"U"),"\n",DimsString(XPre,"X"));
    )
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();
    const bool conjugate = ( orientation == ADJOINT );

    DistMatrixReadProxy<T,T,MC,MR> UProx( UPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<T,MR,  STAR> U12Trans_MR_STAR(g);
    DistMatrix<T,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<T,MC,  STAR> D1_MC_STAR(g);
    
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto U11 = U( IR(k,k+nb), IR(k,k+nb) );
        auto U12 = U( IR(k,k+nb), IR(k+nb,n) );

        auto X1 = X( ALL, IR(k,k+nb) );
        auto X2 = X( ALL, IR(k+nb,n) );

        X1_VC_STAR = X1;
        U11_STAR_STAR = U11;
        LocalTrmm
        ( RIGHT, UPPER, orientation, diag, T(1), U11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
 
        U12Trans_MR_STAR.AlignWith( X2 );
        Transpose( U12, U12Trans_MR_STAR, conjugate );
        D1_MC_STAR.AlignWith( X1 );
        LocalGemm( NORMAL, NORMAL, T(1), X2, U12Trans_MR_STAR, D1_MC_STAR );
        AxpyContract( T(1), D1_MC_STAR, X1 );
    }
}

// Right Upper Adjoint/Transpose (Non)Unit Trmm
//   X := X triu(U)^T, 
//   X := X triu(U)^H,
//   X := X triuu(U)^T, or
//   X := X triuu(U)^H
template<typename T>
void RUT
( Orientation orientation,
  UnitOrNonUnit diag,
  const AbstractDistMatrix<T>& U,
        AbstractDistMatrix<T>& X )
{
    EL_DEBUG_CSE
    // TODO: Come up with a better routing mechanism
    if( U.Height() > 5*X.Height() )
        RUTA( orientation, diag, U, X );
    else
        RUTC( orientation, diag, U, X );
}

} // namespace trmm
} // namespace El

#endif // ifndef EL_TRMM_RUT_HPP
