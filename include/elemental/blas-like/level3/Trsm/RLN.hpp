/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   Copyright (c) 2013, The University of Texas at Austin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TRSM_RLN_HPP
#define ELEM_TRSM_RLN_HPP

#include ELEM_GEMM_INC

namespace elem {
namespace trsm {

// Right Lower Normal (Non)Unit Trsm
//   X := X tril(L)^-1, and
//   X := X trilu(L)^-1
template<typename F>
inline void
RLN
( UnitOrNonUnit diag, const DistMatrix<F>& L, DistMatrix<F>& X,
  bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("trsm::RLN"))
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<F,MR,  STAR> L10Trans_MR_STAR(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,STAR,MC  > X1Trans_STAR_MC(g);
    DistMatrix<F,VC,  STAR> X1_VC_STAR(g);

    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto L10 = LockedViewRange( L, k, 0, k+nb, k    );
        auto L11 = LockedViewRange( L, k, k, k+nb, k+nb );

        auto X0 = ViewRange( X, 0, 0, m, k    );
        auto X1 = ViewRange( X, 0, k, m, k+nb );

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
}

} // namespace trsm
} // namespace elem

#endif // ifndef ELEM_TRSM_RLN_HPP
