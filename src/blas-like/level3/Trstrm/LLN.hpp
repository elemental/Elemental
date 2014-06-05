/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_TRSTRM_LLN_HPP
#define EL_TRSTRM_LLN_HPP

namespace El {
namespace trstrm {

template<typename F>
inline void
LLNUnb( UnitOrNonUnit diag, F alpha, const Matrix<F>& L, Matrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("trstrm::LLNUnb"))
    const bool isUnit = ( diag==UNIT );
    const Int n = L.Height();
    const Int LLDim = L.LDim();
    const Int XLDim = X.LDim();
    const F* LBuffer = L.LockedBuffer();
    F* XBuffer = X.Buffer();

    // X := alpha X
    if( alpha != F(1) )
        for( Int j=0; j<n; ++j ) 
            for( Int i=j; i<n; ++i )
                XBuffer[i+j*XLDim] *= alpha;

    for( Int i=0; i<n; ++i )
    {
        if( !isUnit )
        {
            const F lambda11 = LBuffer[i+i*LLDim];
            for( Int j=0; j<i; ++j )
                XBuffer[i+j*XLDim] /= lambda11;
            XBuffer[i+i*XLDim] /= lambda11;
        }

        const Int l21Height = n - (i+1);
        const F* l21 = &LBuffer[(i+1)+i*LLDim];
        const F* x1L = &XBuffer[i];
        F* X2L = &XBuffer[i+1];
        blas::Geru( l21Height, i+1, F(-1), l21, 1, x1L, XLDim, X2L, XLDim );
    }
}

template<typename F>
inline void
LLN
( UnitOrNonUnit diag, F alpha, const Matrix<F>& L, Matrix<F>& X,
  bool checkIfSingular=true )
{
    DEBUG_ONLY(CallStackEntry cse("trstrm::LLN"))
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
    ScaleTrapezoid( alpha, LOWER, X );
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
        trstrm::LLNUnb( diag, F(1), L11, X11 );
        Gemm( NORMAL, NORMAL, F(-1), L21, X10, F(1), X20 );
        Z11 = X11;
        MakeTriangular( LOWER, Z11 );
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
}

template<typename F>
inline void
LLN
( UnitOrNonUnit diag, F alpha, const DistMatrix<F>& L, DistMatrix<F>& X,
  bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("trstrm::LLN"))
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
    ScaleTrapezoid( alpha, LOWER, X );
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
        LocalTrstrm
        ( LEFT, LOWER, NORMAL, diag, F(1), L11_STAR_STAR, X11_STAR_STAR,
          checkIfSingular );
        X11 = X11_STAR_STAR;
        X11_STAR_MR = X11_STAR_STAR;
        MakeTriangular( LOWER, X11_STAR_MR );

        X10_STAR_MR = X10_STAR_VR;
        X10 = X10_STAR_MR;
        L21_MC_STAR = L21;
        
        LocalGemm
        ( NORMAL, NORMAL, F(-1), L21_MC_STAR, X10_STAR_MR, F(1), X20 );
        LocalGemm
        ( NORMAL, NORMAL, F(-1), L21_MC_STAR, X11_STAR_MR, F(1), X21 );
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
}

} // namespace trstrm
} // namespace El

#endif // ifndef EL_TRSTRM_LLN_HPP
