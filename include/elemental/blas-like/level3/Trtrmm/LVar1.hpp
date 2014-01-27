/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TRTRMM_LVAR1_HPP
#define ELEM_TRTRMM_LVAR1_HPP

#include ELEM_TRMM_INC

namespace elem {
namespace internal {

template<typename T>
inline void
TrtrmmLVar1( Matrix<T>& L, bool conjugate=false )
{
    DEBUG_ONLY(CallStackEntry cse("internal::TrtrmmLVar1"))
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
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
        TrtrmmLUnblocked( L11, conjugate );
        //--------------------------------------------------------------------/

        SlidePartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );
    }
}

template<typename T>
inline void
TrtrmmLVar1( DistMatrix<T>& L, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrtrmmLVar1");
        if( L.Height() != L.Width() )
            LogicError("L must be square");
    )
    const Grid& g = L.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    // Matrix views
    DistMatrix<T>
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

        LocalTrtrmm( LOWER, L11_STAR_STAR, conjugate );
        L11 = L11_STAR_STAR;
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );
    }
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_TRTRMM_LVAR1_HPP
