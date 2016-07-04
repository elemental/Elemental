/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2013, The University of Texas at Austin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRMM_LLT_HPP
#define EL_TRMM_LLT_HPP

namespace El {
namespace trmm {

template<typename T>
void LocalAccumulateLLT
( Orientation orientation,
  UnitOrNonUnit diag,
  T alpha,
  const DistMatrix<T>& L,
  const DistMatrix<T,MC,STAR>& X,
        DistMatrix<T,MR,STAR>& Z )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( L, X, Z );
      if( L.Height() != L.Width() || L.Height() != X.Height() ||
          L.Height() != Z.Height() )
          LogicError
          ("Nonconformal:\n",DimsString(L,"L"),"\n",DimsString(X,"X"),"\n",
           DimsString(Z,"Z"));
      if( X.ColAlign() != L.ColAlign() || Z.ColAlign() != L.RowAlign() )
          LogicError("Partial matrix distributions are misaligned");
    )
    const Int m = Z.Height();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();
    
    DistMatrix<T> D11(g);

    const Int ratio = Max( g.Height(), g.Width() );
    for( Int k=0; k<m; k+=ratio*bsize )
    {
        const Int nb = Min(ratio*bsize,m-k);

        auto L11 = L( IR(k,k+nb), IR(k,k+nb) );
        auto L21 = L( IR(k+nb,m), IR(k,k+nb) );

        auto X1 = X( IR(k,k+nb), ALL );
        auto X2 = X( IR(k+nb,m), ALL );

        auto Z1 = Z( IR(k,k+nb), ALL );

        D11.AlignWith( L11 );
        D11 = L11;
        MakeTrapezoidal( LOWER, D11 );
        if( diag == UNIT )
            FillDiagonal( D11, T(1) );
        LocalGemm( orientation, NORMAL, alpha, D11, X1, T(1), Z1 );
        LocalGemm( orientation, NORMAL, alpha, L21, X2, T(1), Z1 );
    }
}

template<typename T>
void LLTA
( Orientation orientation,
  UnitOrNonUnit diag,
  const AbstractDistMatrix<T>& LPre,
        AbstractDistMatrix<T>& XPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( LPre, XPre );
      if( orientation == NORMAL )
          LogicError("Expected (Conjugate)Transpose option");
      if( LPre.Height() != LPre.Width() || LPre.Height() != XPre.Height() )
          LogicError
          ("Nonconformal: \n",DimsString(LPre,"L"),"\n",DimsString(XPre,"X"))
    )
    const Int m = XPre.Height();
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    DistMatrixReadProxy<T,T,MC,MR> LProx( LPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> XProx( XPre );
    auto& L = LProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<T,MC,STAR> X1_MC_STAR(g);
    DistMatrix<T,MR,STAR> Z1_MR_STAR(g);
    DistMatrix<T,MR,MC  > Z1_MR_MC(g);

    X1_MC_STAR.AlignWith( L );
    Z1_MR_STAR.AlignWith( L );

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto X1 = X( ALL, IR(k,k+nb) );

        X1_MC_STAR = X1;
        Z1_MR_STAR.Resize( m, nb );
        Zero( Z1_MR_STAR );
        LocalAccumulateLLT
        ( orientation, diag, T(1), L, X1_MC_STAR, Z1_MR_STAR );

        Contract( Z1_MR_STAR, Z1_MR_MC );
        X1 = Z1_MR_MC;
    }
}
   
template<typename T>
void LLTCOld
( Orientation orientation,
  UnitOrNonUnit diag,
  const AbstractDistMatrix<T>& LPre,
        AbstractDistMatrix<T>& XPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( LPre, XPre );
      if( orientation == NORMAL )
          LogicError("Expected (Conjugate)Transpose option");
      if( LPre.Height() != LPre.Width() || LPre.Height() != XPre.Height() )
          LogicError
          ("Nonconformal: \n",DimsString(LPre,"L"),"\n",DimsString(XPre,"X"))
    )
    const Int m = XPre.Height();
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();
    const bool conjugate = ( orientation == ADJOINT );

    DistMatrixReadProxy<T,T,MC,MR> LProx( LPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> XProx( XPre );
    auto& L = LProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> D1Trans_MR_STAR(g);
    DistMatrix<T,MR,  MC  > D1Trans_MR_MC(g);
    DistMatrix<T,MC,  MR  > D1(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto L11 = L( IR(k,k+nb), IR(k,k+nb) );
        auto L21 = L( IR(k+nb,m), IR(k,k+nb) );

        auto X1 = X( IR(k,k+nb), ALL );
        auto X2 = X( IR(k+nb,m), ALL );

        X1_STAR_VR = X1;
        L11_STAR_STAR = L11;
        LocalTrmm
        ( LEFT, LOWER, orientation, diag, T(1), L11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
 
        L21_MC_STAR.AlignWith( X2 );
        L21_MC_STAR = L21;
        D1Trans_MR_STAR.AlignWith( X1 );
        LocalGemm
        ( orientation, NORMAL, T(1), X2, L21_MC_STAR, D1Trans_MR_STAR );
        D1Trans_MR_MC.AlignWith( X1 );
        Contract( D1Trans_MR_STAR, D1Trans_MR_MC );
        D1.AlignWith( X1 );
        D1.Resize( nb, n );
        Zero( D1 );
        Transpose( D1Trans_MR_MC.Matrix(), D1.Matrix(), conjugate );
        Axpy( T(1), D1, X1 );
    }
}

template<typename T>
void LLTC
( Orientation orientation,
  UnitOrNonUnit diag,
  const AbstractDistMatrix<T>& LPre,
        AbstractDistMatrix<T>& XPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( LPre, XPre );
      if( orientation == NORMAL )
          LogicError("Expected (Conjugate)Transpose option");
      if( LPre.Height() != LPre.Width() || LPre.Height() != XPre.Height() )
          LogicError
          ("Nonconformal: \n",DimsString(LPre,"L"),"\n",DimsString(XPre,"X"))
    )
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    DistMatrixReadProxy<T,T,MC,MR> LProx( LPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> XProx( XPre );
    auto& L = LProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> X1Trans_MR_STAR(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto L10 = L( IR(k,k+nb), IR(0,k)    );
        auto L11 = L( IR(k,k+nb), IR(k,k+nb) );

        auto X0 = X( IR(0,k),    ALL );
        auto X1 = X( IR(k,k+nb), ALL );

        L10_STAR_MC.AlignWith( X0 );
        L10_STAR_MC = L10;
        X1Trans_MR_STAR.AlignWith( X0 );
        Transpose( X1, X1Trans_MR_STAR );
        LocalGemm
        ( orientation, TRANSPOSE, 
          T(1), L10_STAR_MC, X1Trans_MR_STAR, T(1), X0 );

        L11_STAR_STAR = L11;
        X1_STAR_VR.AlignWith( X1 );
        Transpose( X1Trans_MR_STAR, X1_STAR_VR );
        LocalTrmm
        ( LEFT, LOWER, orientation, diag, T(1), L11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
    }
}

// Left Lower (Conjugate)Transpose (Non)Unit Trmm
//   X := tril(L)^T,
//   X := tril(L)^H,
//   X := trilu(L)^T, or
//   X := trilu(L)^H
template<typename T>
void LLT
( Orientation orientation,
  UnitOrNonUnit diag,
  const AbstractDistMatrix<T>& L,
        AbstractDistMatrix<T>& X )
{
    DEBUG_CSE
    // TODO: Come up with a better routing mechanism
    if( L.Height() > 5*X.Width() )
        LLTA( orientation, diag, L, X );
    else
        LLTC( orientation, diag, L, X );
}

} // namespace trmm
} // namespace El

#endif // ifndef EL_TRMM_LLT_HPP
