/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_TRIANGULARINVERSE_LVAR3_HPP
#define LAPACK_TRIANGULARINVERSE_LVAR3_HPP

#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"

namespace elem {
namespace triangular_inverse {

template<typename F>
inline void
LVar3Unb( UnitOrNonUnit diag, Matrix<F>& L )
{
#ifndef RELEASE
    PushCallStack("triangular_inverse::LVar3Unb");
    if( L.Height() != L.Width() )
        throw std::logic_error("Nonsquare matrices cannot be triangular");
#endif
    const int n = L.Height();
    const int ldl = L.LDim();
    F* LBuffer = L.Buffer();
    for( int j=0; j<n; ++j )
    {
        const F lambda = ( diag==NON_UNIT ? LBuffer[j+j*ldl] : F(1) );
        for( int k=0; k<j; ++k )
            LBuffer[j+k*ldl] /= -lambda;
        blas::Geru
        ( n-(j+1), j, F(1),
          &LBuffer[(j+1)+j*ldl], 1, &LBuffer[j], ldl, 
          &LBuffer[j+1], ldl );
        if( diag == NON_UNIT )
        {
            for( int k=j+1; k<n; ++k )
                LBuffer[k+j*ldl] /= lambda;
            LBuffer[j+j*ldl] = F(1) / LBuffer[j+j*ldl];
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
LVar3( UnitOrNonUnit diag, Matrix<F>& L )
{
#ifndef RELEASE
    PushCallStack("triangular_inverse::LVar3");
    if( L.Height() != L.Width() )
        throw std::logic_error("Nonsquare matrices cannot be triangular");
#endif
    // Matrix views
    Matrix<F> 
        LTL, LTR,  L00, L01, L02,
        LBL, LBR,  L10, L11, L12,
                   L20, L21, L22;

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

        //--------------------------------------------------------------------//
        Trsm( LEFT, LOWER, NORMAL, diag, F(-1), L11, L10 );
        Gemm( NORMAL, NORMAL, F(1), L21, L10, F(1), L20 );
        Trsm( RIGHT, LOWER, NORMAL, diag, F(1), L11, L21 );
        LVar3Unb( diag, L11 );
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

template<typename F>
inline void
LVar3( UnitOrNonUnit diag, DistMatrix<F>& L )
{
#ifndef RELEASE
    PushCallStack("triangular_inverse::LVar3");
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

} // namespace triangular_inverse
} // namespace elem

#endif // ifndef LAPACK_TRIANGULARINVERSE_LVAR3_HPP
