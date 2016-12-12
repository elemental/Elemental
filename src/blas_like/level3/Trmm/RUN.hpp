/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2013, The University of Texas at Austin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRMM_RUN_HPP
#define EL_TRMM_RUN_HPP

namespace El {
namespace trmm {

template<typename T>
void LocalAccumulateRUN
( Orientation orientation,
  UnitOrNonUnit diag, T alpha,
  const DistMatrix<T,MC,  MR  >& U,
  const DistMatrix<T,STAR,MC  >& X,
        DistMatrix<T,MR,  STAR>& ZTrans )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( U, X, ZTrans );
      if( U.Height() != U.Width() || U.Height() != X.Width() ||
          U.Height() != ZTrans.Height() )
          LogicError
          ("Nonconformal:\n",DimsString(U,"U"),"\n",DimsString(X,"X"),"\n",
           DimsString(ZTrans,"Z'"));
      if( X.RowAlign() != U.ColAlign() || ZTrans.ColAlign() != U.RowAlign() )
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

        auto X0 = X( ALL, IR(0,k)    );
        auto X1 = X( ALL, IR(k,k+nb) );
        
        auto Z1Trans = ZTrans( IR(k,k+nb), ALL );

        D11.AlignWith( U11 );
        D11 = U11;
        MakeTrapezoidal( UPPER, D11 );
        if( diag == UNIT )
            FillDiagonal( D11, T(1) );
        LocalGemm( orientation, orientation, alpha, D11, X1, T(1), Z1Trans );
        LocalGemm( orientation, orientation, alpha, U01, X0, T(1), Z1Trans );
    }
}

template<typename T>
void RUNA
( UnitOrNonUnit diag, 
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

    DistMatrixReadProxy<T,T,MC,MR> UProx( UPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<T,STAR,VC  > X1_STAR_VC(g);
    DistMatrix<T,STAR,MC  > X1_STAR_MC(g);
    DistMatrix<T,MR,  STAR> Z1Trans_MR_STAR(g);
    DistMatrix<T,MR,  MC  > Z1Trans_MR_MC(g);

    X1_STAR_VC.AlignWith( U );
    X1_STAR_MC.AlignWith( U );
    Z1Trans_MR_STAR.AlignWith( U );

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto X1 = X( IR(k,k+nb), ALL );

        X1_STAR_VC = X1;
        X1_STAR_MC = X1_STAR_VC;
        Z1Trans_MR_STAR.Resize( X1.Width(), X1.Height() );
        Zero( Z1Trans_MR_STAR );
        LocalAccumulateRUN
        ( TRANSPOSE, diag, T(1), U, X1_STAR_MC, Z1Trans_MR_STAR );

        Z1Trans_MR_MC.AlignWith( X1 );
        Contract( Z1Trans_MR_STAR, Z1Trans_MR_MC );
        Transpose( Z1Trans_MR_MC.Matrix(), X1.Matrix() );
    }
}

template<typename T>
void RUNCOld
( UnitOrNonUnit diag, 
  const AbstractDistMatrix<T>& UPre,
        AbstractDistMatrix<T>& XPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( UPre, XPre );
      if( UPre.Height() != UPre.Width() || XPre.Width() != UPre.Height() )
          LogicError
          ("Nonconformal:\n",DimsString(UPre,"U"),"\n",DimsString(XPre,"X"));
    )
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();

    DistMatrixReadProxy<T,T,MC,MR> UProx( UPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<T,MR,  STAR> U01_MR_STAR(g);
    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g); 
    DistMatrix<T,VC,  STAR> X1_VC_STAR(g);    
    DistMatrix<T,MC,  STAR> D1_MC_STAR(g);
    
    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto U01 = U( IR(0,k),    IR(k,k+nb) );
        auto U11 = U( IR(k,k+nb), IR(k,k+nb) );

        auto X0 = X( ALL, IR(0,k)    );
        auto X1 = X( ALL, IR(k,k+nb) );

        X1_VC_STAR = X1;
        U11_STAR_STAR = U11;
        LocalTrmm
        ( RIGHT, UPPER, NORMAL, diag, T(1), U11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
 
        U01_MR_STAR.AlignWith( X0 );
        U01_MR_STAR = U01;
        D1_MC_STAR.AlignWith( X1 );
        LocalGemm( NORMAL, NORMAL, T(1), X0, U01_MR_STAR, D1_MC_STAR );
        AxpyContract( T(1), D1_MC_STAR, X1 );
    }
}

template<typename T>
void RUNC
( UnitOrNonUnit diag, 
  const AbstractDistMatrix<T>& UPre,
        AbstractDistMatrix<T>& XPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( UPre, XPre );
      if( UPre.Height() != UPre.Width() || XPre.Width() != UPre.Height() )
          LogicError
          ("Nonconformal:\n",DimsString(UPre,"U"),"\n",DimsString(XPre,"X"));
    )
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();

    DistMatrixReadProxy<T,T,MC,MR> UProx( UPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<T,MR,  STAR> U12Trans_MR_STAR(g);
    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<T,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<T,MC,  STAR> X1_MC_STAR(g);
    
    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto U11 = U( IR(k,k+nb), IR(k,k+nb) );
        auto U12 = U( IR(k,k+nb), IR(k+nb,n) );

        auto X1 = X( ALL, IR(k,k+nb) );
        auto X2 = X( ALL, IR(k+nb,n) );

        X1_MC_STAR.AlignWith( X2 );
        X1_MC_STAR = X1;
        U12Trans_MR_STAR.AlignWith( X2 );
        Transpose( U12, U12Trans_MR_STAR );
        LocalGemm
        ( NORMAL, TRANSPOSE, T(1), X1_MC_STAR, U12Trans_MR_STAR, T(1), X2 );

        U11_STAR_STAR = U11;
        X1_VC_STAR.AlignWith( X1 );
        X1_VC_STAR = X1_MC_STAR;
        LocalTrmm
        ( RIGHT, UPPER, NORMAL, diag, T(1), U11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
    }
}

// Right Upper Normal (Non)Unit Trmm
//   X := X triu(U), and
//   X := X triuu(U)
template<typename T>
void RUN
( UnitOrNonUnit diag,
  const AbstractDistMatrix<T>& U,
        AbstractDistMatrix<T>& X )
{
    EL_DEBUG_CSE
    // TODO: Come up with a better routing mechanism
    if( U.Height() > 5*X.Height() )
        RUNA( diag, U, X );
    else
        RUNC( diag, U, X );
}

} // namespace trmm
} // namespace El

#endif // ifndef EL_TRMM_RUN_HPP
