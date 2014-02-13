/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TRTRMM_UVAR1_HPP
#define ELEM_TRTRMM_UVAR1_HPP

#include ELEM_TRMM_INC

namespace elem {
namespace internal {

template<typename T>
inline void
TrtrmmUVar1( Matrix<T>& U, bool conjugate=false )
{
    DEBUG_ONLY(CallStackEntry cse("internal::TrtrmmUVar1"))
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE ); 
    Matrix<T>
        UTL, UTR,  U00, U01, U02,
        UBL, UBR,  U10, U11, U12,
                   U20, U21, U22;

    PartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    while( UTL.Height() < U.Height() && UTL.Width() < U.Height() )
    {
        RepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        //--------------------------------------------------------------------/
        Trrk( UPPER, NORMAL, orientation, T(1), U01, U01, T(1), U00 );
        Trmm( RIGHT, UPPER, orientation, NON_UNIT, T(1), U11, U01 );
        TrtrmmUUnblocked( U11, conjugate );
        //--------------------------------------------------------------------/

        SlidePartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );
    }
}

template<typename T>
inline void
TrtrmmUVar1( DistMatrix<T>& U, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrtrmmUVar1");
        if( U.Height() != U.Width() )
            LogicError("U must be square");
    )
    const Grid& g = U.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    // Matrix views
    DistMatrix<T>
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    // Temporary distributions
    DistMatrix<T,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<T,VC,  STAR> U01_VC_STAR(g);
    DistMatrix<T,VR,  STAR> U01_VR_STAR(g);
    DistMatrix<T,STAR,MR  > U01Trans_STAR_MR(g);
    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g);

    U01_MC_STAR.AlignWith( U );
    U01_VC_STAR.AlignWith( U );
    U01_VR_STAR.AlignWith( U );
    U01Trans_STAR_MR.AlignWith( U );

    PartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    while( UTL.Height() < U.Height() && UTL.Width() < U.Height() )
    {
        RepartitionDownDiagonal 
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        //--------------------------------------------------------------------//
        U01_MC_STAR = U01;
        U01_VC_STAR = U01_MC_STAR;
        U01_VR_STAR = U01_VC_STAR;
        U01_VR_STAR.TransposePartialColAllGather( U01Trans_STAR_MR, conjugate );
        LocalTrrk( UPPER, T(1), U01_MC_STAR, U01Trans_STAR_MR, T(1), U00 );

        U11_STAR_STAR = U11;
        LocalTrmm
        ( RIGHT, UPPER, orientation, NON_UNIT, 
          T(1), U11_STAR_STAR, U01_VC_STAR );
        U01 = U01_VC_STAR;

        LocalTrtrmm( UPPER, U11_STAR_STAR, conjugate );
        U11 = U11_STAR_STAR;
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );
    }
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_TRTRMM_UVAR1_HPP
