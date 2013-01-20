/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   Copyright (c) 2013, The University of Texas at Austin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_TRMM_RLN_HPP
#define BLAS_TRMM_RLN_HPP

namespace elem {
namespace internal {

template<typename T>
inline void
TrmmRLNA
( UnitOrNonUnit diag, 
  T alpha, const DistMatrix<T>& L,
                 DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmRLNA");
    if( L.Grid() != X.Grid() )
        throw std::logic_error("{L,X} must be distributed over the same grid");
#endif
    const Grid& g = L.Grid();

    DistMatrix<T>
        XT(g),  X0(g),
        XB(g),  X1(g),
                X2(g);

    DistMatrix<T,STAR,VC  > X1_STAR_VC(g);
    DistMatrix<T,STAR,MC  > X1_STAR_MC(g);
    DistMatrix<T,MR,  STAR> Z1Trans_MR_STAR(g);
    DistMatrix<T,MR,  MC  > Z1Trans_MR_MC(g);

    X1_STAR_VC.AlignWith( L );
    X1_STAR_MC.AlignWith( L );
    Z1Trans_MR_STAR.AlignWith( L );

    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XT.Height() < X.Height() )
    {
        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );

        Z1Trans_MR_MC.AlignWith( X1 );
        Zeros( X1.Width(), X1.Height(), Z1Trans_MR_STAR );
        //--------------------------------------------------------------------//
        X1_STAR_VC = X1;
        X1_STAR_MC = X1_STAR_VC;
        LocalTrmmAccumulateRLN
        ( TRANSPOSE, diag, alpha, L, X1_STAR_MC, Z1Trans_MR_STAR );

        Z1Trans_MR_MC.SumScatterFrom( Z1Trans_MR_STAR );
        Transpose( Z1Trans_MR_MC.LocalMatrix(), X1.LocalMatrix() );
        //--------------------------------------------------------------------//
        Z1Trans_MR_MC.FreeAlignments();

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
    
template<typename T>
inline void
TrmmRLNCOld
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& L,
                 DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmRLNCOld");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( L.Height() != L.Width() || X.Width() != L.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal TrmmRLNC: \n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<T> XL(g), XR(g),
                  X0(g), X1(g), X2(g);

    // Temporary distributions
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,MR,  STAR> L21_MR_STAR(g);
    DistMatrix<T,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<T,MC,  STAR> D1_MC_STAR(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionRight( X, XL, XR, 0 );
    while( XR.Width() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );
 
        RepartitionRight
        ( XL, /**/ XR,
          X0, /**/ X1, X2 );

        L21_MR_STAR.AlignWith( X2 );
        D1_MC_STAR.AlignWith( X1 );
        Zeros( X1.Height(), X1.Width(), D1_MC_STAR );
        //--------------------------------------------------------------------//
        X1_VC_STAR = X1;
        L11_STAR_STAR = L11;
        LocalTrmm
        ( RIGHT, LOWER, NORMAL, diag, T(1), L11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
 
        L21_MR_STAR = L21;
        LocalGemm( NORMAL, NORMAL, T(1), X2, L21_MR_STAR, T(0), D1_MC_STAR );
        X1.SumScatterUpdate( T(1), D1_MC_STAR );
        //--------------------------------------------------------------------//
        L21_MR_STAR.FreeAlignments();
        D1_MC_STAR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionRight
        ( XL,     /**/ XR,
          X0, X1, /**/ X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
TrmmRLNC
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& L,
                 DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmRLNC");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( L.Height() != L.Width() || X.Width() != L.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal TrmmRLNC: \n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<T> XL(g), XR(g),
                  X0(g), X1(g), X2(g);

    // Temporary distributions
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,MR,  STAR> L10Trans_MR_STAR(g);
    DistMatrix<T,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<T,MC,  STAR> X1_MC_STAR(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionRight( X, XL, XR, 0 );
    while( XR.Width() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );
 
        RepartitionRight
        ( XL, /**/ XR,
          X0, /**/ X1, X2 );

        X1_MC_STAR.AlignWith( X0 );
        L10Trans_MR_STAR.AlignWith( X0 );
        X1_VC_STAR.AlignWith( X1 );
        //--------------------------------------------------------------------//
        X1_MC_STAR = X1;
        L10Trans_MR_STAR.TransposeFrom( L10 );
        LocalGemm
        ( NORMAL, TRANSPOSE, T(1), X1_MC_STAR, L10Trans_MR_STAR, T(1), X0 );

        L11_STAR_STAR = L11;
        X1_VC_STAR = X1_MC_STAR;
        LocalTrmm
        ( RIGHT, LOWER, NORMAL, diag, T(1), L11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
        //--------------------------------------------------------------------//
        X1_MC_STAR.FreeAlignments();
        L10Trans_MR_STAR.FreeAlignments();
        X1_VC_STAR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionRight
        ( XL,     /**/ XR,
          X0, X1, /**/ X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalTrmmAccumulateRLN
( Orientation orientation, UnitOrNonUnit diag, T alpha,
  const DistMatrix<T,MC,  MR  >& L,
  const DistMatrix<T,STAR,MC  >& X_STAR_MC,
        DistMatrix<T,MR,  STAR>& ZAdjOrTrans_MR_STAR )
{
#ifndef RELEASE
    PushCallStack("internal::LocalTrmmAccumulateRLN");
    if( L.Grid() != X_STAR_MC.Grid() || 
        X_STAR_MC.Grid() != ZAdjOrTrans_MR_STAR.Grid() )
        throw std::logic_error
        ("{L,X,Z} must be distributed over the same grid");
    if( L.Height() != L.Width() ||
        L.Height() != X_STAR_MC.Width() ||
        L.Height() != ZAdjOrTrans_MR_STAR.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrmmAccumulateRLN: \n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X[* ,MC] ~ " << X_STAR_MC.Height() << " x "
                               << X_STAR_MC.Width() << "\n"
            << "  Z^H/T[MR,* ] ~ " << ZAdjOrTrans_MR_STAR.Height() << " x "
                                   << ZAdjOrTrans_MR_STAR.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( X_STAR_MC.RowAlignment() != L.ColAlignment() ||
        ZAdjOrTrans_MR_STAR.ColAlignment() != L.RowAlignment() )
        throw std::logic_error("Partial matrix distributions are misaligned");
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<T> D11(g);

    DistMatrix<T,STAR,MC>
        XL_STAR_MC(g), XR_STAR_MC(g),
        X0_STAR_MC(g), X1_STAR_MC(g), X2_STAR_MC(g);

    DistMatrix<T,MR,STAR>
        ZTAdjOrTrans_MR_STAR(g),  Z0AdjOrTrans_MR_STAR(g),
        ZBAdjOrTrans_MR_STAR(g),  Z1AdjOrTrans_MR_STAR(g),
                                   Z2AdjOrTrans_MR_STAR(g);

    const int ratio = std::max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );

    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    LockedPartitionRight( X_STAR_MC,  XL_STAR_MC, XR_STAR_MC, 0 );
    PartitionDown
    ( ZAdjOrTrans_MR_STAR, ZTAdjOrTrans_MR_STAR,
                            ZBAdjOrTrans_MR_STAR, 0 );
    while( LTL.Height() < L.Height() )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        LockedRepartitionRight
        ( XL_STAR_MC, /**/ XR_STAR_MC,
          X0_STAR_MC, /**/ X1_STAR_MC, X2_STAR_MC );

        RepartitionDown
        ( ZTAdjOrTrans_MR_STAR,  Z0AdjOrTrans_MR_STAR,
         /*********************/ /*********************/
                                  Z1AdjOrTrans_MR_STAR,
          ZBAdjOrTrans_MR_STAR,  Z2AdjOrTrans_MR_STAR );

        D11.AlignWith( L11 );
        //--------------------------------------------------------------------//
        D11 = L11;
        MakeTrapezoidal( LEFT, LOWER, 0, D11 );
        if( diag == UNIT )
            SetDiagonalToOne( D11 );
        LocalGemm
        ( orientation, orientation,
          alpha, D11, X1_STAR_MC, T(1), Z1AdjOrTrans_MR_STAR );

        LocalGemm
        ( orientation, orientation,
          alpha, L21, X2_STAR_MC, T(1), Z1AdjOrTrans_MR_STAR );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlideLockedPartitionRight
        ( XL_STAR_MC,             /**/ XR_STAR_MC,
          X0_STAR_MC, X1_STAR_MC, /**/ X2_STAR_MC );

        SlidePartitionDown
        ( ZTAdjOrTrans_MR_STAR,  Z0AdjOrTrans_MR_STAR,
                                 Z1AdjOrTrans_MR_STAR,
         /********************/ /********************/
          ZBAdjOrTrans_MR_STAR,  Z2AdjOrTrans_MR_STAR );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

// Right Lower Normal (Non)Unit Trmm
//   X := X tril(L), and
//   X := X trilu(L)
template<typename T>
inline void
TrmmRLN
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& L,
                 DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmRLN");
#endif
    // TODO: Come up with a better routing mechanism
    if( L.Height() > 5*X.Height() )
        TrmmRLNA( diag, alpha, L, X );
    else
        TrmmRLNC( diag, alpha, L, X );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem

#endif // ifndef BLAS_TRMM_RLN_HPP
