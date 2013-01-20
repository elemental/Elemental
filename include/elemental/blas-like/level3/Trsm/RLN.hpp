/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   Copyright (c) 2013, The University of Texas at Austin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_TRSM_RLN_HPP
#define BLAS_TRSM_RLN_HPP

namespace elem {
namespace internal {

// Right Lower Normal (Non)Unit Trsm
//   X := X tril(L)^-1, and
//   X := X trilu(L)^-1
template<typename F>
inline void
TrsmRLN
( UnitOrNonUnit diag, F alpha, const DistMatrix<F>& L, DistMatrix<F>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("internal::TrsmRLN");
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<F> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<F> XL(g), XR(g),
                  X0(g), X1(g), X2(g);

    // Temporary distributions
    DistMatrix<F,MR,  STAR> L10Trans_MR_STAR(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,STAR,MC  > X1Trans_STAR_MC(g);
    DistMatrix<F,VC,  STAR> X1_VC_STAR(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionUpDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionLeft( X, XL, XR, 0 );
    while( XL.Width() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );
 
        RepartitionLeft
        ( XL,     /**/ XR,
          X0, X1, /**/ X2 );

        X1Trans_STAR_MC.AlignWith( X0 );
        L10Trans_MR_STAR.AlignWith( X0 );
        //--------------------------------------------------------------------//
        L11_STAR_STAR = L11;
        X1_VC_STAR = X1;
        LocalTrsm
        ( RIGHT, LOWER, NORMAL, diag, F(1), L11_STAR_STAR, X1_VC_STAR,
          checkIfSingular );

        // X0[MC,MR] -= X1[MC,* ]   L10[*,MR]
        //            = X1^T[* ,MC] L10^T[MR,* ]
        X1Trans_STAR_MC.TransposeFrom( X1_VC_STAR );
        X1.TransposeFrom( X1Trans_STAR_MC );
        L10Trans_MR_STAR.TransposeFrom( L10 );
        LocalGemm
        ( TRANSPOSE, TRANSPOSE, 
          F(-1), X1Trans_STAR_MC, L10Trans_MR_STAR, F(1), X0 );
        //--------------------------------------------------------------------//
        X1Trans_STAR_MC.FreeAlignments();
        L10Trans_MR_STAR.FreeAlignments();

        SlideLockedPartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        SlidePartitionLeft
        ( XL, /**/ XR,
          X0, /**/ X1, X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem

#endif // ifndef BLAS_TRSM_RLN_HPP
