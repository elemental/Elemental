/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace trsm {

// Left Upper Normal (Non)Unit Trsm
//   X := triu(U)^-1  X, or
//   X := triuu(U)^-1 X

template<typename F>
void LUNLarge
( UnitOrNonUnit diag,
  const AbstractDistMatrix<F>& UPre,
        AbstractDistMatrix<F>& XPre, 
  bool checkIfSingular )
{
    EL_DEBUG_CSE
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<F,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    const Int kLast = LastOffset( m, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        auto X0 = X( ind0, ALL );
        auto X1 = X( ind1, ALL );

        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        X1_STAR_VR    = X1;  // X1[* ,VR] <- X1[MC,MR]
        
        // X1[* ,VR] := U11^-1[* ,* ] X1[* ,VR]
        LocalTrsm
        ( LEFT, UPPER, NORMAL, diag, F(1), U11_STAR_STAR, X1_STAR_VR,
          checkIfSingular );

        X1_STAR_MR.AlignWith( X0 );
        X1_STAR_MR  = X1_STAR_VR; // X1[* ,MR]  <- X1[* ,VR]
        X1          = X1_STAR_MR; // X1[MC,MR] <- X1[* ,MR]
        U01_MC_STAR.AlignWith( X0 );
        U01_MC_STAR = U01;        // U01[MC,* ] <- U01[MC,MR]

        // X0[MC,MR] -= U01[MC,* ] X1[* ,MR]
        LocalGemm( NORMAL, NORMAL, F(-1), U01_MC_STAR, X1_STAR_MR, F(1), X0 );
    }
}

template<typename F>
void LUNMedium
( UnitOrNonUnit diag, 
  const AbstractDistMatrix<F>& UPre,
        AbstractDistMatrix<F>& XPre,
  bool checkIfSingular )
{
    EL_DEBUG_CSE
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<F,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,MR,  STAR> X1Trans_MR_STAR(g);

    const Int kLast = LastOffset( m, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        auto X0 = X( ind0, ALL );
        auto X1 = X( ind1, ALL );

        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        X1Trans_MR_STAR.AlignWith( X0 );
        Transpose( X1, X1Trans_MR_STAR );
        
        // X1^T[MR,* ] := X1^T[MR,* ] U11^-T[* ,* ]
        //              = (U11^-1[* ,* ] X1[* ,MR])^T
        LocalTrsm
        ( RIGHT, UPPER, TRANSPOSE, diag, 
          F(1), U11_STAR_STAR, X1Trans_MR_STAR, checkIfSingular );
        Transpose( X1Trans_MR_STAR, X1 );

        U01_MC_STAR.AlignWith( X0 );
        U01_MC_STAR = U01;  // U01[MC,* ] <- U01[MC,MR]

        // X0[MC,MR] -= U01[MC,* ] X1[* ,MR]
        LocalGemm
        ( NORMAL, TRANSPOSE, F(-1), U01_MC_STAR, X1Trans_MR_STAR, F(1), X0 );
    }
}

template<typename F,Dist colDist>
void LUNSmall
( UnitOrNonUnit diag,
  const DistMatrix<F,colDist,STAR>& U,
        DistMatrix<F,colDist,STAR>& X,
  bool checkIfSingular )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( U, X );
      if( U.Height() != U.Width() || U.Width() != X.Height() )
          LogicError
          ("Nonconformal: \n",DimsString(U,"U"),"\n",DimsString(X,"X"));
      if( U.ColAlign() != X.ColAlign() )
          LogicError("U and X are assumed to be aligned");
    )
    const Int m = X.Height();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g), X1_STAR_STAR(g);

    const Int kLast = LastOffset( m, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        auto X0 = X( ind0, ALL );
        auto X1 = X( ind1, ALL );

        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[VC,* ]
        X1_STAR_STAR = X1;   // X1[* ,* ] <- X1[VC,* ]
        
        // X1[* ,* ] := U11^-1[* ,* ] X1[* ,* ]
        LocalTrsm
        ( LEFT, UPPER, NORMAL, diag,
          F(1), U11_STAR_STAR, X1_STAR_STAR, checkIfSingular );
        X1 = X1_STAR_STAR;

        // X0[VC,* ] -= U01[VC,* ] X1[* ,* ]
        LocalGemm( NORMAL, NORMAL, F(-1), U01, X1_STAR_STAR, F(1), X0 );
    }
}

} // namespace trsm
} // namespace El
