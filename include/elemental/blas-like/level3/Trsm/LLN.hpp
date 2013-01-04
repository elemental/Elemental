/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

// Left Lower NORMAL (Non)Unit Trsm 
//   X := tril(L)^-1  X, or
//   X := trilu(L)^-1 X

// For large numbers of RHS's, e.g., width(X) >> p
template<typename F>
inline void
TrsmLLNLarge
( UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& L, DistMatrix<F>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("internal::TrsmLLNLarge");
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<F> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);
    DistMatrix<F> XT(g),  X0(g), 
                  XB(g),  X1(g),
                          X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );

        L21_MC_STAR.AlignWith( X2 );
        X1_STAR_MR.AlignWith( X2 );
        //--------------------------------------------------------------------//
        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[MC,MR]
        X1_STAR_VR    = X1;  // X1[* ,VR] <- X1[MC,MR]

        // X1[* ,VR] := L11^-1[* ,* ] X1[* ,VR]
        LocalTrsm
        ( LEFT, LOWER, NORMAL, diag, F(1), L11_STAR_STAR, X1_STAR_VR,
          checkIfSingular );

        X1_STAR_MR  = X1_STAR_VR; // X1[* ,MR]  <- X1[* ,VR]
        X1          = X1_STAR_MR; // X1[MC,MR] <- X1[* ,MR]
        L21_MC_STAR = L21;        // L21[MC,* ] <- L21[MC,MR]
        
        // X2[MC,MR] -= L21[MC,* ] X1[* ,MR]
        LocalGemm
        ( NORMAL, NORMAL, F(-1), L21_MC_STAR, X1_STAR_MR, F(1), X2 );
        //--------------------------------------------------------------------//
        L21_MC_STAR.FreeAlignments();
        X1_STAR_MR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12, 
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionDown
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// For medium numbers of RHS's, e.g., width(X) ~= p
template<typename F>
inline void
TrsmLLNMedium
( UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& L, DistMatrix<F>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("internal::TrsmLLNMedium");
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<F> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);
    DistMatrix<F> XT(g),  X0(g), 
                  XB(g),  X1(g),
                          X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<F,MR,  STAR> X1Trans_MR_STAR(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );

        L21_MC_STAR.AlignWith( X2 );
        X1Trans_MR_STAR.AlignWith( X2 );
        //--------------------------------------------------------------------//
        L11_STAR_STAR = L11;                 // L11[* ,* ] <- L11[MC,MR]
        X1Trans_MR_STAR.TransposeFrom( X1 ); // X1[* ,MR] <- X1[MC,MR]

        // X1[* ,MR] := L11^-1[* ,* ] X1[* ,MR]
        // X1^T[MR,* ] := X1^T[MR,* ] L11^-T[* ,* ]
        LocalTrsm
        ( RIGHT, LOWER, TRANSPOSE, diag, 
          F(1), L11_STAR_STAR, X1Trans_MR_STAR, checkIfSingular );

        X1.TransposeFrom( X1Trans_MR_STAR ); // X1[MC,MR]  <- X1[* ,MR]
        L21_MC_STAR = L21;                   // L21[MC,* ] <- L21[MC,MR]
        
        // X2[MC,MR] -= L21[MC,* ] X1[* ,MR]
        LocalGemm
        ( NORMAL, TRANSPOSE, F(-1), L21_MC_STAR, X1Trans_MR_STAR, F(1), X2 );
        //--------------------------------------------------------------------//
        L21_MC_STAR.FreeAlignments();
        X1Trans_MR_STAR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12, 
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionDown
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// For small numbers of RHS's, e.g., width(X) < p
template<typename F>
inline void
TrsmLLNSmall
( UnitOrNonUnit diag,
  F alpha, const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("internal::TrsmLLNSmall");
    if( L.ColAlignment() != X.ColAlignment() )
        throw std::logic_error("L and X are assumed to be aligned");
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<F,VC,STAR> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);
    DistMatrix<F,VC,STAR> XT(g),  X0(g), 
                          XB(g),  X1(g),
                                  X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> X1_STAR_STAR(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );

        //--------------------------------------------------------------------//
        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[VC,* ]
        X1_STAR_STAR = X1;   // X1[* ,* ] <- X1[VC,* ]

        // X1[* ,* ] := (L11[* ,* ])^-1 X1[* ,* ]
        LocalTrsm
        ( LEFT, LOWER, NORMAL, diag, 
          F(1), L11_STAR_STAR, X1_STAR_STAR, checkIfSingular );

        // X2[VC,* ] -= L21[VC,* ] X1[* ,* ]
        LocalGemm( NORMAL, NORMAL, F(-1), L21, X1_STAR_STAR, F(1), X2 );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12, 
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionDown
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
