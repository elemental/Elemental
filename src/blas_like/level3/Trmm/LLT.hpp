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
#ifndef EL_TRMM_LLT_HPP
#define EL_TRMM_LLT_HPP

namespace El {
namespace trmm {

template<typename Ring>
inline void
LocalAccumulateLLT
( Orientation orientation, UnitOrNonUnit diag, Ring alpha,
  const DistMatrix<Ring>& L,
  const DistMatrix<Ring,MC,STAR>& X,
        DistMatrix<Ring,MR,STAR>& Z )
{
    DEBUG_ONLY(
      CSE cse("trmm::LocalAccumulateLLT");
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
    
    DistMatrix<Ring> D11(g);

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
            FillDiagonal( D11, Ring(1) );
        LocalGemm( alpha, D11.Orient(orientation), X1.N(), Ring(1), Z1 );
        LocalGemm( alpha, L21.Orient(orientation), X2.N(), Ring(1), Z1 );
    }
}

template<typename Ring>
inline void
LLTA
( Orientation orientation, UnitOrNonUnit diag,
  const AbstractDistMatrix<Ring>& LPre, AbstractDistMatrix<Ring>& XPre )
{
    DEBUG_ONLY(
      CSE cse("trmm::LLTA");
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

    auto LPtr = ReadProxy<Ring,MC,MR>( &LPre );      auto& L = *LPtr;
    auto XPtr = ReadWriteProxy<Ring,MC,MR>( &XPre ); auto& X = *XPtr;

    DistMatrix<Ring,MC,STAR> X1_MC_STAR(g);
    DistMatrix<Ring,MR,STAR> Z1_MR_STAR(g);
    DistMatrix<Ring,MR,MC  > Z1_MR_MC(g);

    X1_MC_STAR.AlignWith( L );
    Z1_MR_STAR.AlignWith( L );

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto X1 = X( ALL, IR(k,k+nb) );

        X1_MC_STAR = X1;
        Zeros( Z1_MR_STAR, m, nb );
        LocalAccumulateLLT
        ( orientation, diag, Ring(1), L, X1_MC_STAR, Z1_MR_STAR );

        Contract( Z1_MR_STAR, Z1_MR_MC );
        X1 = Z1_MR_MC;
    }
}
   
template<typename Ring>
inline void
LLTCOld
( Orientation orientation, UnitOrNonUnit diag,
  const AbstractDistMatrix<Ring>& LPre, AbstractDistMatrix<Ring>& XPre )
{
    DEBUG_ONLY(
      CSE cse("trmm::LLTCOld");
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

    auto LPtr = ReadProxy<Ring,MC,MR>( &LPre );      auto& L = *LPtr;
    auto XPtr = ReadWriteProxy<Ring,MC,MR>( &XPre ); auto& X = *XPtr;

    DistMatrix<Ring,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<Ring,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<Ring,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<Ring,MR,  STAR> D1Trans_MR_STAR(g);
    DistMatrix<Ring,MR,  MC  > D1Trans_MR_MC(g);
    DistMatrix<Ring,MC,  MR  > D1(g);

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
        ( LEFT, LOWER, orientation, diag, Ring(1), L11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
 
        L21_MC_STAR.AlignWith( X2 );
        L21_MC_STAR = L21;
        D1Trans_MR_STAR.AlignWith( X1 );
        LocalGemm
        ( Ring(1), X2.Orient(orientation), L21_MC_STAR.N(), D1Trans_MR_STAR );
        D1Trans_MR_MC.AlignWith( X1 );
        Contract( D1Trans_MR_STAR, D1Trans_MR_MC );
        D1.AlignWith( X1 );
        Zeros( D1, nb, n );
        Transpose( D1Trans_MR_MC.Matrix(), D1.Matrix(), conjugate );
        X1 += D1;
    }
}

template<typename Ring>
inline void
LLTC
( Orientation orientation, UnitOrNonUnit diag,
  const AbstractDistMatrix<Ring>& LPre, AbstractDistMatrix<Ring>& XPre )
{
    DEBUG_ONLY(
      CSE cse("trmm::LLTC");
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

    auto LPtr = ReadProxy<Ring,MC,MR>( &LPre );      auto& L = *LPtr;
    auto XPtr = ReadWriteProxy<Ring,MC,MR>( &XPre ); auto& X = *XPtr;

    DistMatrix<Ring,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<Ring,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<Ring,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<Ring,MR,  STAR> X1Trans_MR_STAR(g);

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
        ( Ring(1), L10_STAR_MC.Orient(orientation), 
                   X1Trans_MR_STAR.T(), 
          Ring(1), X0 );

        L11_STAR_STAR = L11;
        X1_STAR_VR.AlignWith( X1 );
        Transpose( X1Trans_MR_STAR, X1_STAR_VR );
        LocalTrmm
        ( LEFT, LOWER, orientation, diag, Ring(1), L11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
    }
}

// Left Lower (Conjugate)Transpose (Non)Unit Trmm
//   X := tril(L)^T,
//   X := tril(L)^H,
//   X := trilu(L)^T, or
//   X := trilu(L)^H
template<typename Ring>
inline void
LLT
( Orientation orientation, UnitOrNonUnit diag,
  const AbstractDistMatrix<Ring>& L, AbstractDistMatrix<Ring>& X )
{
    DEBUG_ONLY(CSE cse("trmm::LLT"))
    // TODO: Come up with a better routing mechanism
    if( L.Height() > 5*X.Width() )
        LLTA( orientation, diag, L, X );
    else
        LLTC( orientation, diag, L, X );
}

} // namespace trmm
} // namespace El

#endif // ifndef EL_TRMM_LLT_HPP
