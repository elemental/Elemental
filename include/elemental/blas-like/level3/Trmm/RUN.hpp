/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   Copyright (c) 2013, The University of Texas at Austin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

template<typename T>
inline void
TrmmRUNA
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& U,
                 DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmRUNA");
    if( U.Grid() != X.Grid() )
        throw std::logic_error("{U,X} must be distributed over the same grid");
#endif
    const Grid& g = U.Grid();

    DistMatrix<T>
        XT(g),  X0(g),
        XB(g),  X1(g),
                X2(g);

    DistMatrix<T,STAR,VC  > X1_STAR_VC(g);
    DistMatrix<T,STAR,MC  > X1_STAR_MC(g);
    DistMatrix<T,MR,  STAR> Z1Trans_MR_STAR(g);
    DistMatrix<T,MR,  MC  > Z1Trans_MR_MC(g);

    X1_STAR_VC.AlignWith( U );
    X1_STAR_MC.AlignWith( U );
    Z1Trans_MR_STAR.AlignWith( U );

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
        LocalTrmmAccumulateRUN
        ( TRANSPOSE, diag, alpha, U, X1_STAR_MC, Z1Trans_MR_STAR );

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
TrmmRUNCOld
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& U,
                 DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmRUNCOld");
    if( U.Grid() != X.Grid() )
        throw std::logic_error
        ("U and X must be distributed over the same grid");
    if( U.Height() != U.Width() || X.Width() != U.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal TrmmRUNC: \n"
            << "  U ~ " << U.Height() << " x " << U.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<T> 
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    DistMatrix<T> XL(g), XR(g),
                  X0(g), X1(g), X2(g);

    // Temporary distributions
    DistMatrix<T,MR,  STAR> U01_MR_STAR(g);
    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g); 
    DistMatrix<T,VC,  STAR> X1_VC_STAR(g);    
    DistMatrix<T,MC,  STAR> D1_MC_STAR(g);
    
    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionUpDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    PartitionLeft( X, XL, XR, 0 );
    while( XL.Width() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12, 
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

        RepartitionLeft
        ( XL,     /**/ XR,
          X0, X1, /**/ X2 ); 

        U01_MR_STAR.AlignWith( X0 );
        D1_MC_STAR.AlignWith( X1 );
        Zeros( X1.Height(), X1.Width(), D1_MC_STAR );
        //--------------------------------------------------------------------//
        X1_VC_STAR = X1;
        U11_STAR_STAR = U11;
        LocalTrmm
        ( RIGHT, UPPER, NORMAL, diag, T(1), U11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
 
        U01_MR_STAR = U01;
        LocalGemm( NORMAL, NORMAL, T(1), X0, U01_MR_STAR, T(0), D1_MC_STAR );
        X1.SumScatterUpdate( T(1), D1_MC_STAR );
        //--------------------------------------------------------------------//
        U01_MR_STAR.FreeAlignments();
        D1_MC_STAR.FreeAlignments();

        SlideLockedPartitionUpDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        SlidePartitionLeft
        ( XL, /**/ XR,
          X0, /**/ X1, X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
TrmmRUNC
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& U,
                 DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmRUNC");
    if( U.Grid() != X.Grid() )
        throw std::logic_error
        ("U and X must be distributed over the same grid");
    if( U.Height() != U.Width() || X.Width() != U.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal TrmmRUNC: \n"
            << "  U ~ " << U.Height() << " x " << U.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<T> 
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    DistMatrix<T> XL(g), XR(g),
                  X0(g), X1(g), X2(g);

    // Temporary distributions
    DistMatrix<T,MR,  STAR> U12Trans_MR_STAR(g);
    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<T,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<T,MC,  STAR> X1_MC_STAR(g);
    
    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionUpDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    PartitionLeft( X, XL, XR, 0 );
    while( XL.Width() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12, 
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

        RepartitionLeft
        ( XL,     /**/ XR,
          X0, X1, /**/ X2 ); 

        X1_MC_STAR.AlignWith( X2 );
        U12Trans_MR_STAR.AlignWith( X2 );
        X1_VC_STAR.AlignWith( X1 );
        //--------------------------------------------------------------------//
        X1_MC_STAR = X1;
        U12Trans_MR_STAR.TransposeFrom( U12 );
        LocalGemm
        ( NORMAL, TRANSPOSE, T(1), X1_MC_STAR, U12Trans_MR_STAR, T(1), X2 );

        U11_STAR_STAR = U11;
        X1_VC_STAR = X1_MC_STAR;
        LocalTrmm
        ( RIGHT, UPPER, NORMAL, diag, T(1), U11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
        //--------------------------------------------------------------------//
        X1_MC_STAR.FreeAlignments();
        U12Trans_MR_STAR.FreeAlignments();
        X1_VC_STAR.FreeAlignments();

        SlideLockedPartitionUpDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        SlidePartitionLeft
        ( XL, /**/ XR,
          X0, /**/ X1, X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalTrmmAccumulateRUN
( Orientation orientation, UnitOrNonUnit diag, T alpha,
  const DistMatrix<T,MC,  MR  >& U,
  const DistMatrix<T,STAR,MC  >& X_STAR_MC,
        DistMatrix<T,MR,  STAR>& ZHermOrTrans_MR_STAR )
{
#ifndef RELEASE
    PushCallStack("internal::LocalTrmmAccumulateRUN");
    if( U.Grid() != X_STAR_MC.Grid() ||
        X_STAR_MC.Grid() != ZHermOrTrans_MR_STAR.Grid() )
        throw std::logic_error
        ("{U,X,Z} must be distributed over the same grid");
    if( U.Height() != U.Width() ||
        U.Height() != X_STAR_MC.Width() ||
        U.Height() != ZHermOrTrans_MR_STAR.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrmmAccumulateRUN: \n"
            << "  U ~ " << U.Height() << " x " << U.Width() << "\n"
            << "  X[* ,MC] ~ " << X_STAR_MC.Height() << " x "
                               << X_STAR_MC.Width() << "\n"
            << "  Z^H/T[MR,* ] ~ " << ZHermOrTrans_MR_STAR.Height() << " x "
                                   << ZHermOrTrans_MR_STAR.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( X_STAR_MC.RowAlignment() != U.ColAlignment() ||
        ZHermOrTrans_MR_STAR.ColAlignment() != U.RowAlignment() )
        throw std::logic_error("Partial matrix distributions are misaligned");
#endif
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<T>
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    DistMatrix<T> D11(g);

    DistMatrix<T,STAR,MC>
        XL_STAR_MC(g), XR_STAR_MC(g),
        X0_STAR_MC(g), X1_STAR_MC(g), X2_STAR_MC(g);

    DistMatrix<T,MR,STAR>
        ZTHermOrTrans_MR_STAR(g),  Z0HermOrTrans_MR_STAR(g),
        ZBHermOrTrans_MR_STAR(g),  Z1HermOrTrans_MR_STAR(g),
                                   Z2HermOrTrans_MR_STAR(g);

    const int ratio = std::max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );

    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    LockedPartitionRight( X_STAR_MC,  XL_STAR_MC, XR_STAR_MC, 0 );
    PartitionDown
    ( ZHermOrTrans_MR_STAR, ZTHermOrTrans_MR_STAR,
                            ZBHermOrTrans_MR_STAR, 0 );
    while( UTL.Height() < U.Height() )
    {
        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        LockedRepartitionRight
        ( XL_STAR_MC, /**/ XR_STAR_MC,
          X0_STAR_MC, /**/ X1_STAR_MC, X2_STAR_MC );

        RepartitionDown
        ( ZTHermOrTrans_MR_STAR,  Z0HermOrTrans_MR_STAR,
         /*********************/ /*********************/
                                  Z1HermOrTrans_MR_STAR,
          ZBHermOrTrans_MR_STAR,  Z2HermOrTrans_MR_STAR );

        D11.AlignWith( U11 );
        //--------------------------------------------------------------------//
        D11 = U11;
        MakeTrapezoidal( LEFT, UPPER, 0, D11 );
        if( diag == UNIT )
            SetDiagonalToOne( D11 );
        LocalGemm
        ( orientation, orientation,
          alpha, D11, X1_STAR_MC, T(1), Z1HermOrTrans_MR_STAR );

        LocalGemm
        ( orientation, orientation,
          alpha, U01, X0_STAR_MC, T(1), Z1HermOrTrans_MR_STAR );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

        SlideLockedPartitionRight
        ( XL_STAR_MC,             /**/ XR_STAR_MC,
          X0_STAR_MC, X1_STAR_MC, /**/ X2_STAR_MC );

        SlidePartitionDown
        ( ZTHermOrTrans_MR_STAR,  Z0HermOrTrans_MR_STAR,
                                  Z1HermOrTrans_MR_STAR,
         /*********************/ /*********************/
          ZBHermOrTrans_MR_STAR,  Z2HermOrTrans_MR_STAR );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

// Right Upper Normal (Non)Unit Trmm
//   X := X triu(U), and
//   X := X triuu(U)
template<typename T>
inline void
TrmmRUN
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& U,
                 DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmRUN");
#endif
    // TODO: Come up with a better routing mechanism
    if( U.Height() > 5*X.Height() )
        TrmmRUNA( diag, alpha, U, X );
    else
        TrmmRUNC( diag, alpha, U, X );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
