/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace trsm {

// Right Lower (Conjugate)Transpose (Non)Unit Trsm
//   X := X tril(L)^-T,
//   X := X tril(L)^-H,
//   X := X trilu(L)^-T, or
//   X := X trilu(L)^-H
template<typename F>
inline void
RLT
( Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<F>& L, DistMatrix<F>& X, bool checkIfSingular )
{
    DEBUG_ONLY(
        CallStackEntry cse("trsm::RLT");
        if( orientation == NORMAL )
            LogicError("Expected (Conjugate)Transpose option");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,VR,  STAR> L21_VR_STAR(g);
    DistMatrix<F,STAR,MR  > L21Trans_STAR_MR(g);
    DistMatrix<F,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<F,STAR,MC  > X1Trans_STAR_MC(g);

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto L11 = LockedViewRange( L, k,    k, k+nb, k+nb );
        auto L21 = LockedViewRange( L, k+nb, k, n,    k+nb ); 

        auto X1 = ViewRange( X, 0, k,    m, k+nb );
        auto X2 = ViewRange( X, 0, k+nb, m, n    );

        L11_STAR_STAR = L11; 
        X1_VC_STAR.AlignWith( X2 );
        X1_VC_STAR = X1;  
        
        LocalTrsm
        ( RIGHT, LOWER, orientation, diag, 
          F(1), L11_STAR_STAR, X1_VC_STAR, checkIfSingular );

        X1Trans_STAR_MC.AlignWith( X2 );
        X1_VC_STAR.TransposePartialColAllGather( X1Trans_STAR_MC );
        X1.TransposeRowFilterFrom( X1Trans_STAR_MC );
        L21_VR_STAR.AlignWith( X2 );
        L21_VR_STAR = L21;
        L21Trans_STAR_MR.AlignWith( X2 );
        L21_VR_STAR.TransposePartialColAllGather
        ( L21Trans_STAR_MR, (orientation==ADJOINT) ); 

        // X2[MC,MR] -= X1[MC,*] (L21[MR,*])^(T/H)
        //            = X1^T[* ,MC] (L21^(T/H))[*,MR]
        LocalGemm
        ( TRANSPOSE, NORMAL, 
          F(-1), X1Trans_STAR_MC, L21Trans_STAR_MR, F(1), X2 );
    }
}

} // namespace trsm
} // namespace El
