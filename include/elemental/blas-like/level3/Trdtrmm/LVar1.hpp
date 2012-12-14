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
TrdtrmmLVar1( Orientation orientation, Matrix<F>& L )
{
#ifndef RELEASE
    PushCallStack("internal::TrdtrmmLVar1");
    if( L.Height() != L.Width() )
        throw std::logic_error("L must be square");
    if( orientation == NORMAL )
        throw std::logic_error("Orientation must be (conjugate-)transpose");
#endif
    Matrix<F>
        LTL, LTR,  L00, L01, L02,
        LBL, LBR,  L10, L11, L12,
                   L20, L21, L22;
    Matrix<F> d1, S10;

    PartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    while( LTL.Height() < L.Height() && LTL.Width() < L.Height() )
    {
        RepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        //--------------------------------------------------------------------/
        L11.GetDiagonal( d1 );
        S10 = L10;
        DiagonalSolve( LEFT, NORMAL, d1, L10, true );
        Trrk( LOWER, orientation, NORMAL, F(1), S10, L10, F(1), L00 );
        Trmm( LEFT, LOWER, orientation, UNIT, F(1), L11, L10 );
        TrdtrmmLUnblocked( orientation, L11 );
        //--------------------------------------------------------------------/

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

template<typename F>
inline void
TrdtrmmLVar1( Orientation orientation, DistMatrix<F,MC,MR>& L )
{
#ifndef RELEASE
    PushCallStack("internal::TrdtrmmLVar1");
    if( L.Height() != L.Width() )
        throw std::logic_error("L must be square");
    if( orientation == NORMAL )
        throw std::logic_error("Orientation must be (conjugate-)transpose");
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<F,MC,MR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);
    DistMatrix<F,MD,STAR> d1(g);

    // Temporary distributions
    DistMatrix<F,STAR,VR  > L10_STAR_VR(g);
    DistMatrix<F,STAR,VC  > S10_STAR_VC(g);
    DistMatrix<F,STAR,MC  > S10_STAR_MC(g);
    DistMatrix<F,STAR,MR  > L10_STAR_MR(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);

    L10_STAR_VR.AlignWith( L );
    S10_STAR_VC.AlignWith( L );
    S10_STAR_MC.AlignWith( L );
    L10_STAR_MR.AlignWith( L );

    PartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    while( LTL.Height() < L.Height() && LTL.Width() < L.Height() )
    {
        RepartitionDownDiagonal 
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        //--------------------------------------------------------------------//
        L11.GetDiagonal( d1 );

        L10_STAR_VR = L10;
        S10_STAR_VC = L10_STAR_VR;
        S10_STAR_MC = S10_STAR_VC;
        DiagonalSolve( LEFT, NORMAL, d1, L10_STAR_VR, true );
        L10_STAR_MR = L10_STAR_VR;
        LocalTrrk
        ( LOWER, orientation, F(1), S10_STAR_MC, L10_STAR_MR, F(1), L00 );

        L11_STAR_STAR = L11;
        LocalTrmm
        ( LEFT, LOWER, orientation, UNIT, F(1), L11_STAR_STAR, L10_STAR_VR );
        L10 = L10_STAR_VR;

        LocalTrdtrmm( orientation, LOWER, L11_STAR_STAR );
        L11 = L11_STAR_STAR;
        //--------------------------------------------------------------------//
        d1.FreeAlignments();

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
