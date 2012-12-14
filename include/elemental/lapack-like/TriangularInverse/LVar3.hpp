/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

template<typename F>
inline void
TriangularInverseLVar3( UnitOrNonUnit diag, DistMatrix<F>& L )
{
#ifndef RELEASE
    PushCallStack("internal::TriangularInverseLVar3");
    if( L.Height() != L.Width() )
        throw std::logic_error("Nonsquare matrices cannot be triangular");
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<F> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    // Temporary distributions
    DistMatrix<F,STAR,MR  > L10_STAR_MR(g);
    DistMatrix<F,STAR,VR  > L10_STAR_VR(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<F,VC,  STAR> L21_VC_STAR(g);

    // Start the algorithm
    PartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    while( LTL.Height() < L.Height() )
    {
        RepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        L10_STAR_MR.AlignWith( L20 );
        L21_MC_STAR.AlignWith( L20 );
        //--------------------------------------------------------------------//
        L10_STAR_VR = L10;
        L11_STAR_STAR = L11;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, diag, F(-1), L11_STAR_STAR, L10_STAR_VR );

        L21_MC_STAR = L21;
        L10_STAR_MR = L10_STAR_VR;
        LocalGemm
        ( NORMAL, NORMAL, F(1), L21_MC_STAR, L10_STAR_MR, F(1), L20 );
        L10 = L10_STAR_MR;

        L21_VC_STAR = L21_MC_STAR;
        LocalTrsm
        ( RIGHT, LOWER, NORMAL, diag, F(1), L11_STAR_STAR, L21_VC_STAR );
        LocalTriangularInverse( LOWER, diag, L11_STAR_STAR );
        L11 = L11_STAR_STAR;
        L21 = L21_VC_STAR;
        //--------------------------------------------------------------------//
        L10_STAR_MR.FreeAlignments();
        L21_MC_STAR.FreeAlignments();

        SlidePartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
