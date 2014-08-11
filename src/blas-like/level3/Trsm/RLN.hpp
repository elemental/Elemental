/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   Copyright (c) 2013, The University of Texas at Austin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace trsm {

// Right Lower Normal (Non)Unit Trsm
//   X := X tril(L)^-1, and
//   X := X trilu(L)^-1
template<typename F>
inline void
RLN
( UnitOrNonUnit diag, 
  const AbstractDistMatrix<F>& LPre, AbstractDistMatrix<F>& XPre,
  bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("trsm::RLN"))
    const Int m = XPre.Height();
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    DistMatrix<F> L(g), X(g);
    Copy( LPre, L, READ_PROXY );
    Copy( XPre, X, READ_WRITE_PROXY );

    DistMatrix<F,MR,  STAR> L10Trans_MR_STAR(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,STAR,MC  > X1Trans_STAR_MC(g);
    DistMatrix<F,VC,  STAR> X1_VC_STAR(g);

    const IndexRange outerInd( 0, m );

    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const IndexRange ind0( 0, k    );
        const IndexRange ind1( k, k+nb );

        auto L10 = LockedView( L, ind1, ind0 );
        auto L11 = LockedView( L, ind1, ind1 );

        auto X0 = View( X, outerInd, ind0 );
        auto X1 = View( X, outerInd, ind1 );

        L11_STAR_STAR = L11;
        X1_VC_STAR = X1;
        LocalTrsm
        ( RIGHT, LOWER, NORMAL, diag, F(1), L11_STAR_STAR, X1_VC_STAR,
          checkIfSingular );

        // X0[MC,MR] -= X1[MC,* ]   L10[*,MR]
        //            = X1^T[* ,MC] L10^T[MR,* ]
        X1Trans_STAR_MC.AlignWith( X0 );
        X1_VC_STAR.TransposePartialColAllGather( X1Trans_STAR_MC );
        X1.TransposeRowFilterFrom( X1Trans_STAR_MC );
        L10Trans_MR_STAR.AlignWith( X0 );
        L10.TransposeColAllGather( L10Trans_MR_STAR );
        LocalGemm
        ( TRANSPOSE, TRANSPOSE, 
          F(-1), X1Trans_STAR_MC, L10Trans_MR_STAR, F(1), X0 );
    }
    Copy( X, XPre, RESTORE_READ_WRITE_PROXY );
}

} // namespace trsm
} // namespace El
