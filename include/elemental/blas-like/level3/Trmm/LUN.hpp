/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

namespace elem {
namespace internal {

template<typename T>
inline void
TrmmLUNA
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& U,
                 DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmULNA");
    if( U.Grid() != X.Grid() )
        throw std::logic_error
        ("U and X must be distributed over the same grid");
    if( U.Height() != U.Width() || U.Width() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal TrmmLUNA: \n"
            << "  U ~ " << U.Height() << " x " << U.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = U.Grid();

    DistMatrix<T>
        XL(g), XR(g),
        X0(g), X1(g), X2(g);

    DistMatrix<T,VR,  STAR> X1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > X1Trans_STAR_MR(g);
    DistMatrix<T,MC,  STAR> Z1_MC_STAR(g);

    X1_VR_STAR.AlignWith( U );
    X1Trans_STAR_MR.AlignWith( U );
    Z1_MC_STAR.AlignWith( U );

    PartitionRight
    ( X, XL, XR, 0 );
    while( XL.Width() < X.Width() )
    {
        RepartitionRight
        ( XL, /**/ XR,
          X0, /**/ X1, X2 );

        Zeros( X1.Height(), X1.Width(), Z1_MC_STAR );
        //--------------------------------------------------------------------//
        X1_VR_STAR = X1;
        X1Trans_STAR_MR.TransposeFrom( X1_VR_STAR );
        LocalTrmmAccumulateLUN
        ( TRANSPOSE, diag, alpha, U, X1Trans_STAR_MR, Z1_MC_STAR );

        X1.SumScatterFrom( Z1_MC_STAR );
        //--------------------------------------------------------------------//

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
TrmmLUNC
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& U,
                 DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmLUNC");
    if( U.Grid() != X.Grid() )
        throw std::logic_error
        ("U and X must be distributed over the same grid");
    if( U.Height() != U.Width() || U.Width() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal TrmmLUN: \n"
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
    DistMatrix<T> XT(g),  X0(g),
                  XB(g),  X1(g),
                          X2(g);

    // Temporary distributions
    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<T,STAR,MC  > U12_STAR_MC(g);
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> D1Trans_MR_STAR(g);
    DistMatrix<T,MR,  MC  > D1Trans_MR_MC(g);
    DistMatrix<T,MC,  MR  > D1(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,   U00, /**/ U01, U02,
         /*************/  /******************/
               /**/        U10, /**/ U11, U12,
          UBL, /**/ UBR,   U20, /**/ U21, U22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );

        U12_STAR_MC.AlignWith( X2 );
        D1Trans_MR_STAR.AlignWith( X1 );
        D1Trans_MR_MC.AlignWith( X1 );
        D1.AlignWith( X1 );
        Zeros( X1.Width(), X1.Height(), D1Trans_MR_STAR );
        Zeros( X1.Height(), X1.Width(), D1 );
        //--------------------------------------------------------------------//
        X1_STAR_VR = X1;
        U11_STAR_STAR = U11;
        LocalTrmm
        ( LEFT, UPPER, NORMAL, diag, (T)1, U11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
 
        U12_STAR_MC = U12;
        LocalGemm
        ( TRANSPOSE, TRANSPOSE, (T)1, X2, U12_STAR_MC, (T)0, D1Trans_MR_STAR );
        D1Trans_MR_MC.SumScatterFrom( D1Trans_MR_STAR );
        Transpose( D1Trans_MR_MC.LocalMatrix(), D1.LocalMatrix() );
        Axpy( (T)1, D1, X1 );
       //--------------------------------------------------------------------//
        D1.FreeAlignments();
        D1Trans_MR_MC.FreeAlignments();
        D1Trans_MR_STAR.FreeAlignments();
        U12_STAR_MC.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

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
LocalTrmmAccumulateLUN
( Orientation orientation, UnitOrNonUnit diag, T alpha,
  const DistMatrix<T,MC,  MR  >& U,
  const DistMatrix<T,STAR,MR  >& XAdjOrTrans_STAR_MR,
        DistMatrix<T,MC,  STAR>& Z_MC_STAR )
{
#ifndef RELEASE
    PushCallStack("internal::LocalTrmmAccumulateLUN");
    if( U.Grid() != XAdjOrTrans_STAR_MR.Grid() ||
        XAdjOrTrans_STAR_MR.Grid() != Z_MC_STAR.Grid() )
        throw std::logic_error
        ("{U,X,Z} must be distributed over the same grid");
    if( U.Height() != U.Width() ||
        U.Height() != XAdjOrTrans_STAR_MR.Width() ||
        U.Height() != Z_MC_STAR.Height() ||
        XAdjOrTrans_STAR_MR.Height() != Z_MC_STAR.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrmmAccumulateLUN: \n"
            << "  U ~ " << U.Height() << " x " << U.Width() << "\n"
            << "  X^H/T[* ,MR] ~ " << XAdjOrTrans_STAR_MR.Height() << " x "
                                   << XAdjOrTrans_STAR_MR.Width() << "\n"
            << "  Z[MC,* ] ~ " << Z_MC_STAR.Height() << " x "
                               << Z_MC_STAR.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( XAdjOrTrans_STAR_MR.RowAlignment() != U.RowAlignment() ||
        Z_MC_STAR.ColAlignment() != U.ColAlignment() )
        throw std::logic_error("Partial matrix distributions are misaligned");
#endif
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<T>
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    DistMatrix<T> D11(g);

    DistMatrix<T,STAR,MR>
        XLAdjOrTrans_STAR_MR(g), XRAdjOrTrans_STAR_MR(g),
        X0AdjOrTrans_STAR_MR(g), X1AdjOrTrans_STAR_MR(g), 
        X2AdjOrTrans_STAR_MR(g);

    DistMatrix<T,MC,STAR>
        ZT_MC_STAR(g),  Z0_MC_STAR(g),
        ZB_MC_STAR(g),  Z1_MC_STAR(g),
                        Z2_MC_STAR(g);

    const int ratio = std::max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );

    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    LockedPartitionRight
    ( XAdjOrTrans_STAR_MR, XLAdjOrTrans_STAR_MR, XRAdjOrTrans_STAR_MR, 0 );
    PartitionDown
    ( Z_MC_STAR, ZT_MC_STAR,
                 ZB_MC_STAR, 0 );
    while( UTL.Height() < U.Height() )
    {
        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        LockedRepartitionRight
        ( XLAdjOrTrans_STAR_MR, /**/ XRAdjOrTrans_STAR_MR,
          X0AdjOrTrans_STAR_MR, /**/ X1AdjOrTrans_STAR_MR, X2AdjOrTrans_STAR_MR 
        );

        RepartitionDown
        ( ZT_MC_STAR,  Z0_MC_STAR,
         /**********/ /**********/
                       Z1_MC_STAR,
          ZB_MC_STAR,  Z2_MC_STAR );

        D11.AlignWith( U11 );
        //--------------------------------------------------------------------//
        D11 = U11;
        MakeTrapezoidal( LEFT, UPPER, 0, D11 );
        if( diag == UNIT )
            SetDiagonalToOne( D11 );
        LocalGemm
        ( NORMAL, orientation, alpha, D11, X1AdjOrTrans_STAR_MR,
          (T)1, Z1_MC_STAR );

        LocalGemm
        ( NORMAL, orientation, alpha, U01, X1AdjOrTrans_STAR_MR,
          (T)1, Z0_MC_STAR );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

        SlideLockedPartitionRight
        ( XLAdjOrTrans_STAR_MR,                       /**/ XRAdjOrTrans_STAR_MR,
          X0AdjOrTrans_STAR_MR, X1AdjOrTrans_STAR_MR, /**/ X2AdjOrTrans_STAR_MR 
        );

        SlidePartitionDown
        ( ZT_MC_STAR,  Z0_MC_STAR,
                       Z1_MC_STAR,
         /**********/ /**********/
          ZB_MC_STAR,  Z2_MC_STAR );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

// Left Upper Normal (Non)Unit Trmm
//   X := triu(U)  X, or
//   X := triuu(U) X
template<typename T>
inline void
TrmmLUN
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& U,
                 DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmLUN");
#endif
    // TODO: Come up with a better routing mechanism
    if( U.Height() > 5*X.Width() )
        TrmmLUNA( diag, alpha, U, X );
    else
        TrmmLUNC( diag, alpha, U, X );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
