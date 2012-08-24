/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   Copyright (c) 2012, The University of Texas at Austin
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
TrmmLLNA
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& L,
                 DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmLLNA");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( L.Height() != L.Width() || L.Width() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal TrmmLLNA: \n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = L.Grid();

    DistMatrix<T>
        XL(g), XR(g),
        X0(g), X1(g), X2(g);

    DistMatrix<T,VR,  STAR> X1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > X1Trans_STAR_MR(g);
    DistMatrix<T,MC,  STAR> Z1_MC_STAR(g);

    X1_VR_STAR.AlignWith( L );
    X1Trans_STAR_MR.AlignWith( L );
    Z1_MC_STAR.AlignWith( L );

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
        LocalTrmmAccumulateLLN
        ( TRANSPOSE, diag, alpha, L, X1Trans_STAR_MR, Z1_MC_STAR );

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
TrmmLLNCOld
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& L,
                 DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmLLNCOld");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( L.Height() != L.Width() || L.Width() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal TrmmLLNC: \n"
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
    DistMatrix<T> XT(g),  X0(g),
                  XB(g),  X1(g),
                          X2(g);

    // Temporary distributions
    DistMatrix<T,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> D1Trans_MR_STAR(g);
    DistMatrix<T,MR,  MC  > D1Trans_MR_MC(g);
    DistMatrix<T,MC,  MR  > D1(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionUpDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionUp
    ( X, XT,
         XB, 0 );
    while( XT.Height() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        RepartitionUp
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );

        L10_STAR_MC.AlignWith( X0 );
        D1Trans_MR_STAR.AlignWith( X1 );
        D1Trans_MR_MC.AlignWith( X1 );
        D1.AlignWith( X1 );
        Zeros( X1.Width(), X1.Height(), D1Trans_MR_STAR );
        Zeros( X1.Height(), X1.Width(), D1 );
        //--------------------------------------------------------------------//
        L11_STAR_STAR = L11;
        X1_STAR_VR = X1;
        LocalTrmm
        ( LEFT, LOWER, NORMAL, diag, (T)1, L11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;

        L10_STAR_MC = L10;
        LocalGemm
        ( TRANSPOSE, TRANSPOSE, (T)1, X0, L10_STAR_MC, (T)0, D1Trans_MR_STAR );
        D1Trans_MR_MC.SumScatterFrom( D1Trans_MR_STAR );
        Transpose( D1Trans_MR_MC.LocalMatrix(), D1.LocalMatrix() );
        Axpy( (T)1, D1, X1 );
        //--------------------------------------------------------------------//
        D1.FreeAlignments();
        D1Trans_MR_MC.FreeAlignments();
        D1Trans_MR_STAR.FreeAlignments();
        L10_STAR_MC.FreeAlignments();

        SlideLockedPartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12, 
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        SlidePartitionUp
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
TrmmLLNC
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& L,
                 DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmLLNC");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( L.Height() != L.Width() || L.Width() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal TrmmLLNC: \n"
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
    DistMatrix<T> XT(g),  X0(g),
                  XB(g),  X1(g),
                          X2(g);

    // Temporary distributions
    DistMatrix<T,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> X1Trans_MR_STAR(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionUpDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionUp
    ( X, XT,
         XB, 0 );
    while( XT.Height() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        RepartitionUp
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );

        L21_MC_STAR.AlignWith( X2 );
        X1Trans_MR_STAR.AlignWith( X2 );
        X1_STAR_VR.AlignWith( X1 );
        //--------------------------------------------------------------------//
        L21_MC_STAR = L21;
        X1Trans_MR_STAR.TransposeFrom( X1 );
        LocalGemm
        ( NORMAL, TRANSPOSE, (T)1, L21_MC_STAR, X1Trans_MR_STAR, (T)1, X2 );

        L11_STAR_STAR = L11;
        X1_STAR_VR.TransposeFrom( X1Trans_MR_STAR );
        LocalTrmm
        ( LEFT, LOWER, NORMAL, diag, (T)1, L11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
        //--------------------------------------------------------------------//
        L21_MC_STAR.FreeAlignments();
        X1Trans_MR_STAR.FreeAlignments();
        X1_STAR_VR.FreeAlignments();

        SlideLockedPartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12, 
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        SlidePartitionUp
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalTrmmAccumulateLLN
( Orientation orientation, UnitOrNonUnit diag, T alpha,
  const DistMatrix<T,MC,  MR  >& L,
  const DistMatrix<T,STAR,MR  >& XAdjOrTrans_STAR_MR,
        DistMatrix<T,MC,  STAR>& Z_MC_STAR )
{
#ifndef RELEASE
    PushCallStack("internal::LocalTrmmAccumulateLLN");
    if( L.Grid() != XAdjOrTrans_STAR_MR.Grid() ||
        XAdjOrTrans_STAR_MR.Grid() != Z_MC_STAR.Grid() )
        throw std::logic_error
        ("{L,X,Z} must be distributed over the same grid");
    if( L.Height() != L.Width() ||
        L.Height() != XAdjOrTrans_STAR_MR.Width() ||
        L.Height() != Z_MC_STAR.Height() ||
        XAdjOrTrans_STAR_MR.Height() != Z_MC_STAR.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrmmAccumulateLLN: \n" 
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X^H/T[* ,MR] ~ " << XAdjOrTrans_STAR_MR.Height() << " x "
                                   << XAdjOrTrans_STAR_MR.Width() << "\n"
            << "  Z[MC,* ] ~ " << Z_MC_STAR.Height() << " x "
                               << Z_MC_STAR.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( XAdjOrTrans_STAR_MR.RowAlignment() != L.RowAlignment() ||
        Z_MC_STAR.ColAlignment() != L.ColAlignment() )
        throw std::logic_error("Partial matrix distributions are misaligned");
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

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
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    LockedPartitionRight
    ( XAdjOrTrans_STAR_MR, XLAdjOrTrans_STAR_MR, XRAdjOrTrans_STAR_MR, 0 );
    PartitionDown
    ( Z_MC_STAR, ZT_MC_STAR,
                 ZB_MC_STAR, 0 );
    while( LTL.Height() < L.Height() )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        LockedRepartitionRight
        ( XLAdjOrTrans_STAR_MR, /**/ XRAdjOrTrans_STAR_MR,
          X0AdjOrTrans_STAR_MR, /**/ X1AdjOrTrans_STAR_MR, 
                                     X2AdjOrTrans_STAR_MR );

        RepartitionDown
        ( ZT_MC_STAR,  Z0_MC_STAR,
         /**********/ /**********/
                       Z1_MC_STAR,
          ZB_MC_STAR,  Z2_MC_STAR );

        D11.AlignWith( L11 );
        //--------------------------------------------------------------------//
        D11 = L11;
        MakeTrapezoidal( LEFT, LOWER, 0, D11 );
        if( diag == UNIT )
            SetDiagonalToOne( D11 );
        LocalGemm
        ( NORMAL, orientation, alpha, D11, X1AdjOrTrans_STAR_MR,
          (T)1, Z1_MC_STAR );

        LocalGemm
        ( NORMAL, orientation, alpha, L21, X1AdjOrTrans_STAR_MR,
          (T)1, Z2_MC_STAR );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

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

// Left Lower Normal (Non)Unit Trmm 
//   X := tril(L)  X, or
//   X := trilu(L) X
template<typename T>
inline void
TrmmLLN
( UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& L,
                 DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmLLN");
#endif
    // TODO: Come up with a better routing mechanism
    if( L.Height() > 5*X.Width() )
        TrmmLLNA( diag, alpha, L, X );
    else
        TrmmLLNC( diag, alpha, L, X );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
