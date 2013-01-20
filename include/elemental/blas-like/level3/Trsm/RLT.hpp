/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_TRSM_RLT_HPP
#define BLAS_TRSM_RLT_HPP

namespace elem {
namespace internal {

// Right Lower (Conjugate)Transpose (Non)Unit Trsm
//   X := X tril(L)^-T,
//   X := X tril(L)^-H,
//   X := X trilu(L)^-T, or
//   X := X trilu(L)^-H
template<typename F>
inline void
TrsmRLT
( Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& L, DistMatrix<F>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("internal::TrsmRLT");
    if( orientation == NORMAL )
        throw std::logic_error("TrsmRLT expects a (Conjugate)Transpose option");
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
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,VR,  STAR> L21_VR_STAR(g);
    DistMatrix<F,STAR,MR  > L21AdjOrTrans_STAR_MR(g);
    DistMatrix<F,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<F,STAR,MC  > X1Trans_STAR_MC(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionRight( X, XL, XR, 0 );
    while( XR.Width() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        RepartitionRight
        ( XL, /**/     XR,
          X0, /**/ X1, X2 );

        X1_VC_STAR.AlignWith( X2 );
        X1Trans_STAR_MC.AlignWith( X2 );
        L21_VR_STAR.AlignWith( X2 );
        L21AdjOrTrans_STAR_MR.AlignWith( X2 );
        //--------------------------------------------------------------------//
        L11_STAR_STAR = L11; 
        X1_VC_STAR = X1;  
        
        LocalTrsm
        ( RIGHT, LOWER, orientation, diag, 
          F(1), L11_STAR_STAR, X1_VC_STAR, checkIfSingular );

        X1Trans_STAR_MC.TransposeFrom( X1_VC_STAR );
        X1.TransposeFrom( X1Trans_STAR_MC );
        L21_VR_STAR = L21;
        if( orientation == ADJOINT )
            L21AdjOrTrans_STAR_MR.AdjointFrom( L21_VR_STAR );
        else
            L21AdjOrTrans_STAR_MR.TransposeFrom( L21_VR_STAR );

        // X2[MC,MR] -= X1[MC,*] (L21[MR,*])^(T/H)
        //            = X1^T[* ,MC] (L21^(T/H))[*,MR]
        LocalGemm
        ( TRANSPOSE, NORMAL, 
          F(-1), X1Trans_STAR_MC, L21AdjOrTrans_STAR_MR, F(1), X2 );
        //--------------------------------------------------------------------//
        X1_VC_STAR.FreeAlignments();
        X1Trans_STAR_MC.FreeAlignments();
        L21_VR_STAR.FreeAlignments();
        L21AdjOrTrans_STAR_MR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionRight
        ( XL,     /**/ XR,
          X0, X1, /**/ X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem

#endif // ifndef BLAS_TRSM_RLT_HPP
