/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_TRTRSM_LLN_HPP
#define BLAS_TRTRSM_LLN_HPP

namespace elem {
namespace internal {

template<typename F>
inline void
TrtrsmLLNUnb( UnitOrNonUnit diag, F alpha, const Matrix<F>& L, Matrix<F>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrtrsmLLNUnb");
#endif
    const bool isUnit = ( diag==UNIT );
    const int n = L.Height();
    const int LLDim = L.LDim();
    const int XLDim = X.LDim();
    const F* LBuffer = L.LockedBuffer();
    F* XBuffer = X.Buffer();

    // X := alpha X
    if( alpha != F(1) )
        for( int j=0; j<n; ++j ) 
            for( int i=j; i<n; ++i )
                XBuffer[i+j*XLDim] *= alpha;

    for( int i=0; i<n; ++i )
    {
        if( !isUnit )
        {
            const F lambda11 = LBuffer[i+i*LLDim];
            for( int j=0; j<i; ++j )
                XBuffer[i+j*XLDim] /= lambda11;
            XBuffer[i+i*XLDim] /= lambda11;
        }

        const int l21Height = n - (i+1);
        const F* l21 = &LBuffer[(i+1)+i*LLDim];
        const F* x1L = &XBuffer[i];
        F* X2L = &XBuffer[i+1];
        blas::Geru( l21Height, i+1, F(-1), l21, 1, x1L, XLDim, X2L, XLDim );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
TrtrsmLLN
( UnitOrNonUnit diag, F alpha, const Matrix<F>& L, Matrix<F>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("internal::TrtrsmLLN");
#endif
    // Matrix views
    Matrix<F> 
        LTL, LTR,  L00, L01, L02,
        LBL, LBR,  L10, L11, L12,
                   L20, L21, L22;
    Matrix<F> 
        XTL, XTR,  X00, X01, X02,
        XBL, XBR,  X10, X11, X12,
                   X20, X21, X22;

    Matrix<F> Z11;

    // Start the algorithm
    ScaleTrapezoid( alpha, LEFT, LOWER, 0, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDownDiagonal
    ( X, XTL, XTR,
         XBL, XBR, 0 );
    while( XBR.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        RepartitionDownDiagonal
        ( XTL, /**/ XTR,  X00, /**/ X01, X02,
         /*************/ /******************/
               /**/       X10, /**/ X11, X12,
          XBL, /**/ XBR,  X20, /**/ X21, X22 );

        //--------------------------------------------------------------------//
        Trsm( LEFT, LOWER, NORMAL, diag, F(1), L11, X10, checkIfSingular );
        TrtrsmLLNUnb( diag, F(1), L11, X11 );
        Gemm( NORMAL, NORMAL, F(-1), L21, X10, F(1), X20 );
        Z11 = X11;
        MakeTrapezoidal( LEFT, LOWER, 0, Z11 );
        Gemm( NORMAL, NORMAL, F(-1), L21, Z11, F(1), X21 );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12, 
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionDownDiagonal
        ( XTL, /**/ XTR,  X00, X01, /**/ X02,
               /**/       X10, X11, /**/ X12, 
         /*************/ /******************/
          XBL, /**/ XBR,  X20, X21, /**/ X22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
TrtrsmLLN
( UnitOrNonUnit diag, F alpha, const DistMatrix<F>& L, DistMatrix<F>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("internal::TrtrsmLLN");
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<F> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);
    DistMatrix<F> 
        XTL(g), XTR(g),  X00(g), X01(g), X02(g),
        XBL(g), XBR(g),  X10(g), X11(g), X12(g),
                         X20(g), X21(g), X22(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<F,STAR,MR  > X10_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X10_STAR_VR(g);
    DistMatrix<F,STAR,MR  > X11_STAR_MR(g);
    DistMatrix<F,STAR,STAR> X11_STAR_STAR(g);

    // Start the algorithm
    ScaleTrapezoid( alpha, LEFT, LOWER, 0, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDownDiagonal
    ( X, XTL, XTR,
         XBL, XBR, 0 );
    while( XBR.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        RepartitionDownDiagonal
        ( XTL, /**/ XTR,  X00, /**/ X01, X02,
         /*************/ /******************/
               /**/       X10, /**/ X11, X12,
          XBL, /**/ XBR,  X20, /**/ X21, X22 );

        L21_MC_STAR.AlignWith( X20 );
        X10_STAR_MR.AlignWith( X20 );
        X11_STAR_MR.AlignWith( X21 );
        //--------------------------------------------------------------------//
        L11_STAR_STAR = L11; 
        X11_STAR_STAR = X11; 
        X10_STAR_VR = X10;

        LocalTrsm
        ( LEFT, LOWER, NORMAL, diag, F(1), L11_STAR_STAR, X10_STAR_VR,
          checkIfSingular );
        LocalTrtrsm
        ( LEFT, LOWER, NORMAL, diag, F(1), L11_STAR_STAR, X11_STAR_STAR,
          checkIfSingular );
        X11 = X11_STAR_STAR;
        X11_STAR_MR = X11_STAR_STAR;
        MakeTrapezoidal( LEFT, LOWER, 0, X11_STAR_MR );

        X10_STAR_MR = X10_STAR_VR;
        X10 = X10_STAR_MR;
        L21_MC_STAR = L21;
        
        LocalGemm
        ( NORMAL, NORMAL, F(-1), L21_MC_STAR, X10_STAR_MR, F(1), X20 );
        LocalGemm
        ( NORMAL, NORMAL, F(-1), L21_MC_STAR, X11_STAR_MR, F(1), X21 );
        //--------------------------------------------------------------------//
        L21_MC_STAR.FreeAlignments();
        X10_STAR_MR.FreeAlignments();
        X11_STAR_MR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12, 
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionDownDiagonal
        ( XTL, /**/ XTR,  X00, X01, /**/ X02,
               /**/       X10, X11, /**/ X12, 
         /*************/ /******************/
          XBL, /**/ XBR,  X20, X21, /**/ X22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem

#endif // ifndef BLAS_TRTRSM_LLN_HPP
