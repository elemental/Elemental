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
#ifndef EL_TRMM_RUT_HPP
#define EL_TRMM_RUT_HPP

namespace El {
namespace trmm {

template<typename Ring>
inline void
LocalAccumulateRUT
( UnitOrNonUnit diag, Ring alpha,
  const DistMatrix<Ring>& U,
  const DistMatrix<Ring,MR,STAR>& XTrans,
        DistMatrix<Ring,MC,STAR>& ZTrans )
{
    DEBUG_ONLY(
      CSE cse("trmm::LocalAccumulateRUT");
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

    DistMatrix<Ring> D11(g);

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
            FillDiagonal( D11, Ring(1) );
        LocalGemm( alpha, D11.N(), X1Trans.N(), Ring(1), Z1Trans );
        LocalGemm( alpha, U01.N(), X1Trans.N(), Ring(1), Z0Trans );
    }
}

template<typename Ring>
inline void
RUTA
( Orientation orientation, UnitOrNonUnit diag,
  const AbstractDistMatrix<Ring>& UPre, AbstractDistMatrix<Ring>& XPre )
{
    DEBUG_ONLY(
      CSE cse("trmm::RUTA");
      AssertSameGrids( UPre, XPre );
      // TODO: More input checks
    )
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();
    const bool conjugate = ( orientation == ADJOINT );

    auto UPtr = ReadProxy<Ring,MC,MR>( &UPre );      auto& U = *UPtr;
    auto XPtr = ReadWriteProxy<Ring,MC,MR>( &XPre ); auto& X = *XPtr;

    DistMatrix<Ring,MR,  STAR> X1Trans_MR_STAR(g);
    DistMatrix<Ring,MC,  STAR> Z1Trans_MC_STAR(g);
    DistMatrix<Ring,MC,  MR  > Z1Trans(g);
    DistMatrix<Ring,MR,  MC  > Z1Trans_MR_MC(g);

    X1Trans_MR_STAR.AlignWith( U );
    Z1Trans_MC_STAR.AlignWith( U );

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto X1 = X( IR(k,k+nb), ALL );

        Transpose( X1, X1Trans_MR_STAR, conjugate );
        Zeros( Z1Trans_MC_STAR, X.Width(), nb );
        LocalAccumulateRUT
        ( diag, Ring(1), U, X1Trans_MR_STAR, Z1Trans_MC_STAR );

        Contract( Z1Trans_MC_STAR, Z1Trans );
        Z1Trans_MR_MC.AlignWith( X1 );
        Z1Trans_MR_MC = Z1Trans;
        Transpose( Z1Trans_MR_MC.Matrix(), X1.Matrix(), conjugate );
    }
}

template<typename Ring>
inline void
RUTC
( Orientation orientation, UnitOrNonUnit diag,
  const AbstractDistMatrix<Ring>& UPre, AbstractDistMatrix<Ring>& XPre )
{
    DEBUG_ONLY(
      CSE cse("trmm::RUTC");
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

    auto UPtr = ReadProxy<Ring,MC,MR>( &UPre );      auto& U = *UPtr;
    auto XPtr = ReadWriteProxy<Ring,MC,MR>( &XPre ); auto& X = *XPtr;

    DistMatrix<Ring,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<Ring,MR,  STAR> U12Trans_MR_STAR(g);
    DistMatrix<Ring,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<Ring,MC,  STAR> D1_MC_STAR(g);
    
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
        ( RIGHT, UPPER, orientation, diag, Ring(1), U11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
 
        U12Trans_MR_STAR.AlignWith( X2 );
        Transpose( U12, U12Trans_MR_STAR, conjugate );
        D1_MC_STAR.AlignWith( X1 );
        LocalGemm( Ring(1), X2.N(), U12Trans_MR_STAR.N(), D1_MC_STAR );
        AxpyContract( Ring(1), D1_MC_STAR, X1 );
    }
}

// Right Upper Adjoint/Transpose (Non)Unit Trmm
//   X := X triu(U)^T, 
//   X := X triu(U)^H,
//   X := X triuu(U)^T, or
//   X := X triuu(U)^H
template<typename Ring>
inline void
RUT
( Orientation orientation, UnitOrNonUnit diag,
  const AbstractDistMatrix<Ring>& U, AbstractDistMatrix<Ring>& X )
{
    DEBUG_ONLY(CSE cse("trmm::RUT"))
    // TODO: Come up with a better routing mechanism
    if( U.Height() > 5*X.Height() )
        RUTA( orientation, diag, U, X );
    else
        RUTC( orientation, diag, U, X );
}

} // namespace trmm
} // namespace El

#endif // ifndef EL_TRMM_RUT_HPP
