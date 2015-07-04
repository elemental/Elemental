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
#ifndef EL_TRMM_RUN_HPP
#define EL_TRMM_RUN_HPP

namespace El {
namespace trmm {

template<typename Ring>
inline void
LocalAccumulateRUN
( Orientation orient, UnitOrNonUnit diag, Ring alpha,
  const DistMatrix<Ring,MC,  MR  >& U,
  const DistMatrix<Ring,STAR,MC  >& X,
        DistMatrix<Ring,MR,  STAR>& ZTrans )
{
    DEBUG_ONLY(
      CSE cse("trmm::LocalAccumulateRUN");
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

    DistMatrix<Ring> D11(g);

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
            FillDiagonal( D11, Ring(1) );
        LocalGemm
        ( alpha, D11.Orient(orient), X1.Orient(orient), Ring(1), Z1Trans );
        LocalGemm
        ( alpha, U01.Orient(orient), X0.Orient(orient), Ring(1), Z1Trans );
    }
}

template<typename Ring>
inline void
RUNA
( UnitOrNonUnit diag, 
  const AbstractDistMatrix<Ring>& UPre, AbstractDistMatrix<Ring>& XPre )
{
    DEBUG_ONLY(
      CSE cse("trmm::RUNA");
      AssertSameGrids( UPre, XPre );
      // TODO: More input checks
    )
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();

    auto UPtr = ReadProxy<Ring,MC,MR>( &UPre );      auto& U = *UPtr;
    auto XPtr = ReadWriteProxy<Ring,MC,MR>( &XPre ); auto& X = *XPtr;

    DistMatrix<Ring,STAR,VC  > X1_STAR_VC(g);
    DistMatrix<Ring,STAR,MC  > X1_STAR_MC(g);
    DistMatrix<Ring,MR,  STAR> Z1Trans_MR_STAR(g);
    DistMatrix<Ring,MR,  MC  > Z1Trans_MR_MC(g);

    X1_STAR_VC.AlignWith( U );
    X1_STAR_MC.AlignWith( U );
    Z1Trans_MR_STAR.AlignWith( U );

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto X1 = X( IR(k,k+nb), ALL );

        X1_STAR_VC = X1;
        X1_STAR_MC = X1_STAR_VC;
        Zeros( Z1Trans_MR_STAR, X1.Width(), X1.Height() );
        LocalAccumulateRUN
        ( TRANSPOSE, diag, Ring(1), U, X1_STAR_MC, Z1Trans_MR_STAR );

        Z1Trans_MR_MC.AlignWith( X1 );
        Contract( Z1Trans_MR_STAR, Z1Trans_MR_MC );
        Transpose( Z1Trans_MR_MC.Matrix(), X1.Matrix() );
    }
}

template<typename Ring>
inline void
RUNCOld
( UnitOrNonUnit diag, 
  const AbstractDistMatrix<Ring>& UPre, AbstractDistMatrix<Ring>& XPre )
{
    DEBUG_ONLY(
      CSE cse("trmm::RUNCOld");
      AssertSameGrids( UPre, XPre );
      if( UPre.Height() != UPre.Width() || XPre.Width() != UPre.Height() )
          LogicError
          ("Nonconformal:\n",DimsString(UPre,"U"),"\n",DimsString(XPre,"X"));
    )
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();

    auto UPtr = ReadProxy<Ring,MC,MR>( &UPre );      auto& U = *UPtr;
    auto XPtr = ReadWriteProxy<Ring,MC,MR>( &XPre ); auto& X = *XPtr;

    DistMatrix<Ring,MR,  STAR> U01_MR_STAR(g);
    DistMatrix<Ring,STAR,STAR> U11_STAR_STAR(g); 
    DistMatrix<Ring,VC,  STAR> X1_VC_STAR(g);    
    DistMatrix<Ring,MC,  STAR> D1_MC_STAR(g);
    
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
        ( RIGHT, UPPER, NORMAL, diag, Ring(1), U11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
 
        U01_MR_STAR.AlignWith( X0 );
        U01_MR_STAR = U01;
        D1_MC_STAR.AlignWith( X1 );
        LocalGemm( Ring(1), X0.N(), U01_MR_STAR.N(), D1_MC_STAR );
        AxpyContract( Ring(1), D1_MC_STAR, X1 );
    }
}

template<typename Ring>
inline void
RUNC
( UnitOrNonUnit diag, 
  const AbstractDistMatrix<Ring>& UPre, AbstractDistMatrix<Ring>& XPre )
{
    DEBUG_ONLY(
      CSE cse("trmm::RUNC");
      AssertSameGrids( UPre, XPre );
      if( UPre.Height() != UPre.Width() || XPre.Width() != UPre.Height() )
          LogicError
          ("Nonconformal:\n",DimsString(UPre,"U"),"\n",DimsString(XPre,"X"));
    )
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();

    auto UPtr = ReadProxy<Ring,MC,MR>( &UPre );      auto& U = *UPtr;
    auto XPtr = ReadWriteProxy<Ring,MC,MR>( &XPre ); auto& X = *XPtr;

    DistMatrix<Ring,MR,  STAR> U12Trans_MR_STAR(g);
    DistMatrix<Ring,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<Ring,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<Ring,MC,  STAR> X1_MC_STAR(g);
    
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
        LocalGemm( Ring(1), X1_MC_STAR.N(), U12Trans_MR_STAR.T(), Ring(1), X2 );

        U11_STAR_STAR = U11;
        X1_VC_STAR.AlignWith( X1 );
        X1_VC_STAR = X1_MC_STAR;
        LocalTrmm
        ( RIGHT, UPPER, NORMAL, diag, Ring(1), U11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
    }
}

// Right Upper Normal (Non)Unit Trmm
//   X := X triu(U), and
//   X := X triuu(U)
template<typename Ring>
inline void
RUN
( UnitOrNonUnit diag, 
  const AbstractDistMatrix<Ring>& U, AbstractDistMatrix<Ring>& X )
{
    DEBUG_ONLY(CSE cse("trmm::RUN"))
    // TODO: Come up with a better routing mechanism
    if( U.Height() > 5*X.Height() )
        RUNA( diag, U, X );
    else
        RUNC( diag, U, X );
}

} // namespace trmm
} // namespace El

#endif // ifndef EL_TRMM_RUN_HPP
