/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

template<typename T>
inline void
TrtrmmLVar1( Orientation orientation, Matrix<T>& L )
{
#ifndef RELEASE
    PushCallStack("internal::TrtrmmLVar1");
    if( orientation == NORMAL )
        throw std::logic_error("Must be (conjugate-)transposed");
#endif
    Matrix<T>
        LTL, LTR,  L00, L01, L02,
        LBL, LBR,  L10, L11, L12,
                   L20, L21, L22;

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
        Trrk( LOWER, orientation, NORMAL, T(1), L10, L10, T(1), L00 );
        Trmm( LEFT, LOWER, orientation, NON_UNIT, T(1), L11, L10 );
        TrtrmmLUnblocked( orientation, L11 );
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

template<typename T>
inline void
TrtrmmLVar1( Orientation orientation, DistMatrix<T,MC,MR>& L )
{
#ifndef RELEASE
    PushCallStack("internal::TrtrmmLVar1");
    if( L.Height() != L.Width() )
        throw std::logic_error("L must be square");
    if( orientation == NORMAL )
        throw std::logic_error("Must be (conjugate-)transposed");
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T,MC,MR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    // Temporary distributions
    DistMatrix<T,STAR,VR  > L10_STAR_VR(g);
    DistMatrix<T,STAR,VC  > L10_STAR_VC(g);
    DistMatrix<T,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<T,STAR,MR  > L10_STAR_MR(g);
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);

    L10_STAR_VR.AlignWith( L );
    L10_STAR_VC.AlignWith( L );
    L10_STAR_MC.AlignWith( L );
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
        L10_STAR_VR = L10;
        L10_STAR_VC = L10_STAR_VR;
        L10_STAR_MC = L10_STAR_VC;
        L10_STAR_MR = L10_STAR_VR;
        LocalTrrk
        ( LOWER, orientation, T(1), L10_STAR_MC, L10_STAR_MR, T(1), L00 );

        L11_STAR_STAR = L11;
        LocalTrmm
        ( LEFT, LOWER, orientation, NON_UNIT, 
          T(1), L11_STAR_STAR, L10_STAR_VR );
        L10 = L10_STAR_VR;

        LocalTrtrmm( orientation, LOWER, L11_STAR_STAR );
        L11 = L11_STAR_STAR;
        //--------------------------------------------------------------------//

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
