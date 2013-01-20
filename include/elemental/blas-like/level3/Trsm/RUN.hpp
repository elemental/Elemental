/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_TRSM_RUN_HPP
#define BLAS_TRSM_RUN_HPP

namespace elem {
namespace internal {

// Right Upper Normal (Non)Unit Trsm
//   X := X triu(U)^-1, and
//   X := X triuu(U)^-1
template<typename F>
inline void
TrsmRUN
( UnitOrNonUnit diag, F alpha, const DistMatrix<F>& U, DistMatrix<F>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("internal::TrsmRUN");
#endif
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<F> 
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    DistMatrix<F> XL(g), XR(g),
                  X0(g), X1(g), X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g); 
    DistMatrix<F,STAR,MR  > U12_STAR_MR(g);
    DistMatrix<F,VC,  STAR> X1_VC_STAR(g);    
    DistMatrix<F,STAR,MC  > X1Trans_STAR_MC(g);
    
    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    PartitionRight( X, XL, XR, 0 );
    while( XR.Width() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12, 
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        RepartitionRight
        ( XL, /**/     XR,
          X0, /**/ X1, X2 ); 

        X1_VC_STAR.AlignWith( X2 );
        X1Trans_STAR_MC.AlignWith( X2 );
        U12_STAR_MR.AlignWith( X2 );
        //--------------------------------------------------------------------//
        U11_STAR_STAR = U11; 
        X1_VC_STAR = X1;

        LocalTrsm
        ( RIGHT, UPPER, NORMAL, diag, F(1), U11_STAR_STAR, X1_VC_STAR,
          checkIfSingular );

        X1Trans_STAR_MC.TransposeFrom( X1_VC_STAR );
        X1.TransposeFrom( X1Trans_STAR_MC );
        U12_STAR_MR = U12; 

        // X2[MC,MR] -= X1[MC,* ] U12[* ,MR]
        //            = X1^T[* ,MC] U12[* ,MR]
        LocalGemm
        ( TRANSPOSE, NORMAL, F(-1), X1Trans_STAR_MC, U12_STAR_MR, F(1), X2 );
        //--------------------------------------------------------------------//
        X1_VC_STAR.FreeAlignments();
        X1Trans_STAR_MC.FreeAlignments();
        U12_STAR_MR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

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

#endif // ifndef BLAS_TRSM_RUN_HPP
