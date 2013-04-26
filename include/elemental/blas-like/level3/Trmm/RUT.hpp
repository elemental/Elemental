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
#ifndef BLAS_TRMM_RUT_HPP
#define BLAS_TRMM_RUT_HPP

#include "elemental/blas-like/level1/Axpy.hpp"
#include "elemental/blas-like/level1/MakeTriangular.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level1/SetDiagonal.hpp"
#include "elemental/blas-like/level1/Transpose.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {
namespace internal {

template<typename T>
inline void
LocalTrmmAccumulateRUT
( UnitOrNonUnit diag, T alpha,
  const DistMatrix<T>& U,
  const DistMatrix<T,MR,STAR>& XTrans_MR_STAR,
        DistMatrix<T,MC,STAR>& ZTrans_MC_STAR )
{
#ifndef RELEASE
    CallStackEntry entry("internal::LocalTrmmAccumulateRUT");
    if( U.Grid() != XTrans_MR_STAR.Grid() ||
        XTrans_MR_STAR.Grid() != ZTrans_MC_STAR.Grid() )
        throw std::logic_error
        ("{U,X,Z} must be distributed over the same grid");
    if( U.Height() != U.Width() ||
        U.Height() != XTrans_MR_STAR.Height() ||
        U.Height() != ZTrans_MC_STAR.Height() ||
        XTrans_MR_STAR.Width() != ZTrans_MC_STAR.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrmmAccumulateRUT: \n"
            << "  U ~ " << U.Height() << " x " << U.Width() << "\n"
            << "  X^H/T[MR,* ] ~ " << XTrans_MR_STAR.Height() << " x "
                                   << XTrans_MR_STAR.Width() << "\n"
            << "  Z^H/T[MC,* ] ~ " << ZTrans_MC_STAR.Height() << " x "
                                   << ZTrans_MC_STAR.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( XTrans_MR_STAR.ColAlignment() != U.RowAlignment() ||
        ZTrans_MC_STAR.ColAlignment() != U.ColAlignment() )
        throw std::logic_error("Partial matrix distributions are misaligned");
#endif
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<T>
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    DistMatrix<T> D11(g);

    DistMatrix<T,MR,STAR>
        XTTrans_MR_STAR(g),  X0Trans_MR_STAR(g),
        XBTrans_MR_STAR(g),  X1Trans_MR_STAR(g),
                             X2Trans_MR_STAR(g);

    DistMatrix<T,MC,STAR>
        ZTTrans_MC_STAR(g),  Z0Trans_MC_STAR(g),
        ZBTrans_MC_STAR(g),  Z1Trans_MC_STAR(g),
                             Z2Trans_MC_STAR(g);

    const int ratio = std::max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );

    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    LockedPartitionDown
    ( XTrans_MR_STAR, XTTrans_MR_STAR,
                      XBTrans_MR_STAR, 0 );
    PartitionDown
    ( ZTrans_MC_STAR, ZTTrans_MC_STAR,
                      ZBTrans_MC_STAR, 0 );
    while( UTL.Height() < U.Height() )
    {
        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        LockedRepartitionDown
        ( XTTrans_MR_STAR,  X0Trans_MR_STAR,
         /***************/ /***************/
                            X1Trans_MR_STAR,
          XBTrans_MR_STAR,  X2Trans_MR_STAR );

        RepartitionDown
        ( ZTTrans_MC_STAR,  Z0Trans_MC_STAR,
         /***************/ /***************/
                            Z1Trans_MC_STAR,
          ZBTrans_MC_STAR,  Z2Trans_MC_STAR );

        D11.AlignWith( U11 );
        //--------------------------------------------------------------------//
        D11 = U11;
        MakeTriangular( UPPER, D11 );
        if( diag == UNIT )
            SetDiagonal( D11, T(1) );
        LocalGemm
        ( NORMAL, NORMAL, alpha, D11, X1Trans_MR_STAR, T(1), Z1Trans_MC_STAR );
        LocalGemm
        ( NORMAL, NORMAL, alpha, U01, X1Trans_MR_STAR, T(1), Z0Trans_MC_STAR );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

        SlideLockedPartitionDown
        ( XTTrans_MR_STAR,  X0Trans_MR_STAR,
                            X1Trans_MR_STAR,
         /***************/ /***************/
          XBTrans_MR_STAR,  X2Trans_MR_STAR );

        SlidePartitionDown
        ( ZTTrans_MC_STAR,  Z0Trans_MC_STAR,
                            Z1Trans_MC_STAR,
         /***************/ /***************/
          ZBTrans_MC_STAR,  Z2Trans_MC_STAR );
    }
    PopBlocksizeStack();
}

template<typename T>
inline void
TrmmRUTA
( Orientation orientation, UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& U,
                 DistMatrix<T>& X )
{
#ifndef RELEASE
    CallStackEntry entry("internal::TrmmRUTA");
    if( U.Grid() != X.Grid() )
        throw std::logic_error("{U,X} must be distributed over the same grid");
#endif
    const Grid& g = U.Grid();
    const bool conjugate = ( orientation == ADJOINT );

    DistMatrix<T>
        XT(g),  X0(g),
        XB(g),  X1(g),
                X2(g);

    DistMatrix<T,MR,  STAR> X1Trans_MR_STAR(g);
    DistMatrix<T,MC,  STAR> Z1Trans_MC_STAR(g);
    DistMatrix<T,MC,  MR  > Z1Trans(g);
    DistMatrix<T,MR,  MC  > Z1Trans_MR_MC(g);

    X1Trans_MR_STAR.AlignWith( U );
    Z1Trans_MC_STAR.AlignWith( U );

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
        //--------------------------------------------------------------------//
        X1Trans_MR_STAR.TransposeFrom( X1, conjugate );
        Zeros( Z1Trans_MC_STAR, X1.Width(), X1.Height() );
        LocalTrmmAccumulateRUT
        ( diag, alpha,U, X1Trans_MR_STAR, Z1Trans_MC_STAR );

        Z1Trans.SumScatterFrom( Z1Trans_MC_STAR );
        Z1Trans_MR_MC = Z1Trans;
        Transpose( Z1Trans_MR_MC.Matrix(), X1.Matrix(), conjugate );
        //--------------------------------------------------------------------//
        Z1Trans_MR_MC.FreeAlignments();

        SlidePartitionDown
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );
    }
}

template<typename T>
inline void
TrmmRUTC
( Orientation orientation, 
  UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& U,
                 DistMatrix<T>& X )
{
#ifndef RELEASE
    CallStackEntry entry("internal::TrmmRUTC");
    if( U.Grid() != X.Grid() )
        throw std::logic_error
        ("U and X must be distributed over the same grid");
    if( orientation == NORMAL )
        throw std::logic_error("TrmmRUTC expects an Adjoint/Transpose option");
    if( U.Height() != U.Width() || X.Width() != U.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal TrmmRUTC: \n"
            << "  U ~ " << U.Height() << " x " << U.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = U.Grid();
    const bool conjugate = ( orientation == ADJOINT );

    // Matrix views
    DistMatrix<T> 
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    DistMatrix<T> XL(g), XR(g),
                  X0(g), X1(g), X2(g);

    // Temporary distributions
    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<T,MR,  STAR> U12Trans_MR_STAR(g);
    DistMatrix<T,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<T,MC,  STAR> D1_MC_STAR(g);
    
    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    PartitionRight( X, XL, XR, 0 );
    while( XR.Width() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        RepartitionRight
        ( XL, /**/ XR,
          X0, /**/ X1, X2 );

        U12Trans_MR_STAR.AlignWith( X2 );
        D1_MC_STAR.AlignWith( X1 );
        //--------------------------------------------------------------------//
        X1_VC_STAR = X1;
        U11_STAR_STAR = U11;
        LocalTrmm
        ( RIGHT, UPPER, orientation, diag, T(1), U11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
 
        U12Trans_MR_STAR.TransposeFrom( U12, conjugate );
        LocalGemm( NORMAL, NORMAL, T(1), X2, U12Trans_MR_STAR, D1_MC_STAR );
        X1.SumScatterUpdate( T(1), D1_MC_STAR );
        //--------------------------------------------------------------------//
        U12Trans_MR_STAR.FreeAlignments();
        D1_MC_STAR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12, 
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

        SlidePartitionRight
        ( XL,     /**/ XR,
          X0, X1, /**/ X2 );
    }
}

// Right Upper Adjoint/Transpose (Non)Unit Trmm
//   X := X triu(U)^T, 
//   X := X triu(U)^H,
//   X := X triuu(U)^T, or
//   X := X triuu(U)^H
template<typename T>
inline void
TrmmRUT
( Orientation orientation, 
  UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& U,
                 DistMatrix<T>& X )
{
#ifndef RELEASE
    CallStackEntry entry("internal::TrmmRUT");
#endif
    // TODO: Come up with a better routing mechanism
    if( U.Height() > 5*X.Height() )
        TrmmRUTA( orientation, diag, alpha, U, X );
    else
        TrmmRUTC( orientation, diag, alpha, U, X );
}

} // namespace internal
} // namespace elem

#endif // ifndef BLAS_TRMM_RUT_HPP
