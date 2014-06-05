/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace trsm {

// Left Lower NORMAL (Non)Unit Trsm 
//   X := tril(L)^-1  X, or
//   X := trilu(L)^-1 X

// For large numbers of RHS's, e.g., width(X) >> p
template<typename F>
inline void
LLNLarge
( UnitOrNonUnit diag, const DistMatrix<F>& L, DistMatrix<F>& X, 
  bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("trsm::LLNLarge"))
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto L11 = LockedViewRange( L, k,    k, k+nb, k+nb );
        auto L21 = LockedViewRange( L, k+nb, k, m,    k+nb );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[MC,MR]
        X1_STAR_VR    = X1;  // X1[* ,VR] <- X1[MC,MR]

        // X1[* ,VR] := L11^-1[* ,* ] X1[* ,VR]
        LocalTrsm
        ( LEFT, LOWER, NORMAL, diag, F(1), L11_STAR_STAR, X1_STAR_VR,
          checkIfSingular );

        X1_STAR_MR.AlignWith( X2 );
        X1_STAR_MR  = X1_STAR_VR; // X1[* ,MR]  <- X1[* ,VR]
        X1          = X1_STAR_MR; // X1[MC,MR] <- X1[* ,MR]
        L21_MC_STAR.AlignWith( X2 );
        L21_MC_STAR = L21;        // L21[MC,* ] <- L21[MC,MR]
        
        // X2[MC,MR] -= L21[MC,* ] X1[* ,MR]
        LocalGemm( NORMAL, NORMAL, F(-1), L21_MC_STAR, X1_STAR_MR, F(1), X2 );
    }
}

// For medium numbers of RHS's, e.g., width(X) ~= p
template<typename F>
inline void
LLNMedium
( UnitOrNonUnit diag, const DistMatrix<F>& L, DistMatrix<F>& X, 
  bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("trsm::LLNMedium"))
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<F,MR,  STAR> X1Trans_MR_STAR(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto L11 = LockedViewRange( L, k,    k, k+nb, k+nb );
        auto L21 = LockedViewRange( L, k+nb, k, m,    k+nb );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[MC,MR]
        X1Trans_MR_STAR.AlignWith( X2 );
        X1.TransposeColAllGather( X1Trans_MR_STAR ); // X1[* ,MR] <- X1[MC,MR]

        // X1^T[MR,* ] := X1^T[MR,* ] L11^-T[* ,* ]
        //              = (L11^-1[* ,* ] X1[* ,MR])^T
        LocalTrsm
        ( RIGHT, LOWER, TRANSPOSE, diag, 
          F(1), L11_STAR_STAR, X1Trans_MR_STAR, checkIfSingular );

        X1.TransposeColFilterFrom( X1Trans_MR_STAR ); // X1[MC,MR] <- X1[* ,MR]
        L21_MC_STAR.AlignWith( X2 );
        L21_MC_STAR = L21;                   // L21[MC,* ] <- L21[MC,MR]
        
        // X2[MC,MR] -= L21[MC,* ] X1[* ,MR]
        LocalGemm
        ( NORMAL, TRANSPOSE, F(-1), L21_MC_STAR, X1Trans_MR_STAR, F(1), X2 );
    }
}

// For small numbers of RHS's, e.g., width(X) < p
template<typename F,Dist colDist>
inline void
LLNSmall
( UnitOrNonUnit diag, 
  const DistMatrix<F,colDist,STAR>& L, 
        DistMatrix<F,colDist,STAR>& X, bool checkIfSingular )
{
    DEBUG_ONLY(
        CallStackEntry cse("trsm::LLNSmall");
        if( L.ColAlign() != X.ColAlign() )
            LogicError("L and X are assumed to be aligned");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g), X1_STAR_STAR(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto L11 = LockedViewRange( L, k,    k, k+nb, k+nb );
        auto L21 = LockedViewRange( L, k+nb, k, m,    k+nb );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[VC,* ]
        X1_STAR_STAR = X1;   // X1[* ,* ] <- X1[VC,* ]

        // X1[* ,* ] := (L11[* ,* ])^-1 X1[* ,* ]
        LocalTrsm
        ( LEFT, LOWER, NORMAL, diag, 
          F(1), L11_STAR_STAR, X1_STAR_STAR, checkIfSingular );

        // X2[VC,* ] -= L21[VC,* ] X1[* ,* ]
        LocalGemm( NORMAL, NORMAL, F(-1), L21, X1_STAR_STAR, F(1), X2 );
    }
}

} // namespace trsm
} // namespace El
