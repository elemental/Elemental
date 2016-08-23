/*
   Copyright (c) 2009-2016, Jack Poulson
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
void RLT
( Orientation orientation,
  UnitOrNonUnit diag,
  const AbstractDistMatrix<F>& LPre,
        AbstractDistMatrix<F>& XPre, 
  bool checkIfSingular )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( orientation == NORMAL )
          LogicError("Expected (Conjugate)Transpose option");
    )
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> LProx( LPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& L = LProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,VR,  STAR> L21_VR_STAR(g);
    DistMatrix<F,STAR,MR  > L21Trans_STAR_MR(g);
    DistMatrix<F,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<F,STAR,MC  > X1Trans_STAR_MC(g);

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto L11 = L( ind1, ind1 );
        auto L21 = L( ind2, ind1 );

        auto X1 = X( ALL, ind1 );
        auto X2 = X( ALL, ind2 );

        L11_STAR_STAR = L11; 
        X1_VC_STAR.AlignWith( X2 );
        X1_VC_STAR = X1;  
        
        LocalTrsm
        ( RIGHT, LOWER, orientation, diag, 
          F(1), L11_STAR_STAR, X1_VC_STAR, checkIfSingular );

        X1Trans_STAR_MC.AlignWith( X2 );
        Transpose( X1_VC_STAR, X1Trans_STAR_MC );
        Transpose( X1Trans_STAR_MC, X1 );
        L21_VR_STAR.AlignWith( X2 );
        L21_VR_STAR = L21;
        L21Trans_STAR_MR.AlignWith( X2 );
        Transpose( L21_VR_STAR, L21Trans_STAR_MR, (orientation==ADJOINT) );

        // X2[MC,MR] -= X1[MC,*] (L21[MR,*])^(T/H)
        //            = X1^T[* ,MC] (L21^(T/H))[*,MR]
        LocalGemm
        ( TRANSPOSE, NORMAL, 
          F(-1), X1Trans_STAR_MC, L21Trans_STAR_MR, F(1), X2 );
    }
}

} // namespace trsm
} // namespace El
