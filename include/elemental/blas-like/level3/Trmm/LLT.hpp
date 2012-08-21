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
TrmmLLTA
( Orientation orientation, 
  UnitOrNonUnit diag,
  T alpha, 
  const DistMatrix<T>& L,
        DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmLLTA");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( orientation == NORMAL )
        throw std::logic_error
        ("TrmmLLTA expects a (Conjugate)Transpose option");
    if( L.Height() != L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal TrmmLLTA: \n"
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
    DistMatrix<T>
        XL(g), XR(g),
        X0(g), X1(g), X2(g);

    DistMatrix<T,MC,STAR> X1_MC_STAR(g);
    DistMatrix<T,MR,STAR> Z1_MR_STAR(g);
    DistMatrix<T,MR,MC  > Z1_MR_MC(g);

    X1_MC_STAR.AlignWith( L );
    Z1_MR_STAR.AlignWith( L );

    PartitionRight
    ( X, XL, XR, 0 );
    while( XL.Width() < X.Width() )
    {
        RepartitionRight
        ( XL, /**/ XR,
          X0, /**/ X1, X2 );

        Zeros( X1.Height(), X1.Width(), Z1_MR_STAR );
        //--------------------------------------------------------------------//
        X1_MC_STAR = X1;
        LocalTrmmAccumulateLLT
        ( orientation, diag, alpha, L, X1_MC_STAR, Z1_MR_STAR );

        Z1_MR_MC.SumScatterFrom( Z1_MR_STAR );
        X1 = Z1_MR_MC;
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
TrmmLLTCOld
( Orientation orientation, 
  UnitOrNonUnit diag,
  T alpha, 
  const DistMatrix<T>& L,
        DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmLLTCOld");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( orientation == NORMAL )
        throw std::logic_error("TrmmLLT expects a (Conjugate)Transpose option");
    if( L.Height() != L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal TrmmLLTC: \n"
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
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> D1AdjOrTrans_MR_STAR(g);
    DistMatrix<T,MR,  MC  > D1AdjOrTrans_MR_MC(g);
    DistMatrix<T,MC,  MR  > D1(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 ); 

        L21_MC_STAR.AlignWith( X2 );
        D1AdjOrTrans_MR_STAR.AlignWith( X1 );
        D1AdjOrTrans_MR_MC.AlignWith( X1 );
        D1.AlignWith( X1 );
        Zeros( X1.Width(), X1.Height(), D1AdjOrTrans_MR_STAR );
        Zeros( X1.Height(), X1.Width(), D1 );
        //--------------------------------------------------------------------//
        X1_STAR_VR = X1;
        L11_STAR_STAR = L11;
        LocalTrmm
        ( LEFT, LOWER, orientation, diag, (T)1, L11_STAR_STAR, X1_STAR_VR );
        X1 = X1_STAR_VR;
 
        L21_MC_STAR = L21;
        LocalGemm
        ( orientation, NORMAL, 
          (T)1, X2, L21_MC_STAR, (T)0, D1AdjOrTrans_MR_STAR );
        D1AdjOrTrans_MR_MC.SumScatterFrom( D1AdjOrTrans_MR_STAR );
        if( orientation == TRANSPOSE )
            Transpose( D1AdjOrTrans_MR_MC.LocalMatrix(), D1.LocalMatrix() );
        else
            Adjoint( D1AdjOrTrans_MR_MC.LocalMatrix(), D1.LocalMatrix() );
        Axpy( (T)1, D1, X1 );
        //--------------------------------------------------------------------//
        D1.FreeAlignments();
        D1AdjOrTrans_MR_MC.FreeAlignments();
        D1AdjOrTrans_MR_STAR.FreeAlignments();
        L21_MC_STAR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

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
TrmmLLTC
( Orientation orientation, 
  UnitOrNonUnit diag,
  T alpha, 
  const DistMatrix<T>& L,
        DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmLLTC");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( orientation == NORMAL )
        throw std::logic_error("TrmmLLT expects a (Conjugate)Transpose option");
    if( L.Height() != L.Width() || L.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal TrmmLLTC: \n"
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
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<T,STAR,VR  > X1_STAR_VR(g);
    DistMatrix<T,MR,  STAR> X1Trans_MR_STAR(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 ); 

        L10_STAR_MC.AlignWith( X0 );
        X1Trans_MR_STAR.AlignWith( X0 );
        X1_STAR_VR.AlignWith( X1 );
        //--------------------------------------------------------------------//
        L10_STAR_MC = L10;
        X1Trans_MR_STAR.TransposeFrom( X1 );
        LocalGemm
        ( orientation, TRANSPOSE, 
          (T)1, L10_STAR_MC, X1Trans_MR_STAR, (T)1, X0 );

        L11_STAR_STAR = L11;
        X1_STAR_VR.TransposeFrom( X1Trans_MR_STAR );
        LocalTrmm( LEFT, LOWER, orientation, NON_UNIT, (T)1, L11_STAR_STAR );
        X1 = X1_STAR_VR;
        //--------------------------------------------------------------------//
        L10_STAR_MC.FreeAlignments();
        X1Trans_MR_STAR.FreeAlignments();
        X1_STAR_VR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

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
LocalTrmmAccumulateLLT
( Orientation orientation, UnitOrNonUnit diag, T alpha,
  const DistMatrix<T>& L,
  const DistMatrix<T,MC,STAR>& X_MC_STAR,
        DistMatrix<T,MR,STAR>& Z_MR_STAR )
{
#ifndef RELEASE
    PushCallStack("internal::LocalTrmmAccumulateLLT");
    if( L.Grid() != X_MC_STAR.Grid() ||
        X_MC_STAR.Grid() != Z_MR_STAR.Grid() )
        throw std::logic_error
        ("{L,X,Z} must be distributed over the same grid");
    if( L.Height() != L.Width() ||
        L.Height() != X_MC_STAR.Height() ||
        L.Height() != Z_MR_STAR.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrmmAccumulateLLT: " << "\n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X[MC,* ] ~ " << X_MC_STAR.Height() << " x "
                               << X_MC_STAR.Width() << "\n"
            << "  Z[MR,* ] ` " << Z_MR_STAR.Height() << " x "
                               << Z_MR_STAR.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( X_MC_STAR.ColAlignment() != L.ColAlignment() ||
        Z_MR_STAR.ColAlignment() != L.RowAlignment() )
        throw std::logic_error("Partial matrix distributions are misaligned");
#endif
    const Grid& g = L.Grid();
    
    // Matrix views
    DistMatrix<T>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<T> D11(g);

    DistMatrix<T,MC,STAR>
        XT_MC_STAR(g),  X0_MC_STAR(g),
        XB_MC_STAR(g),  X1_MC_STAR(g),
                        X2_MC_STAR(g);

    DistMatrix<T,MR,STAR>
        ZT_MR_STAR(g),  Z0_MR_STAR(g),
        ZB_MR_STAR(g),  Z1_MR_STAR(g),
                        Z2_MR_STAR(g);

    const int ratio = std::max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );

    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    LockedPartitionDown
    ( X_MC_STAR, XT_MC_STAR,
                 XB_MC_STAR, 0 );
    PartitionDown
    ( Z_MR_STAR, ZT_MR_STAR,
                 ZB_MR_STAR, 0 );
    while( LTL.Height() < L.Height() )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        LockedRepartitionDown
        ( XT_MC_STAR,  X0_MC_STAR,
         /**********/ /**********/
                       X1_MC_STAR,
          XB_MC_STAR,  X2_MC_STAR );

        RepartitionDown
        ( ZT_MR_STAR,  Z0_MR_STAR,
         /**********/ /**********/
                       Z1_MR_STAR,
          ZB_MR_STAR,  Z2_MR_STAR );

        D11.AlignWith( L11 );
        //--------------------------------------------------------------------//
        D11 = L11;
        MakeTrapezoidal( LEFT, LOWER, 0, D11 );
        if( diag == UNIT )
            SetDiagonalToOne( D11 );
        LocalGemm
        ( orientation, NORMAL,
          alpha, D11, X1_MC_STAR, (T)1, Z1_MR_STAR );

        LocalGemm
        ( orientation, NORMAL,
          alpha, L21, X2_MC_STAR, (T)1, Z1_MR_STAR );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlideLockedPartitionDown
        ( XT_MC_STAR,  X0_MC_STAR,
                       X1_MC_STAR,
         /**********/ /**********/
          XB_MC_STAR,  X2_MC_STAR );

        SlidePartitionDown
        ( ZT_MR_STAR,  Z0_MR_STAR,
                       Z1_MR_STAR,
         /**********/ /**********/
          ZB_MR_STAR,  Z2_MR_STAR );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

// Left Lower (Conjugate)Transpose (Non)Unit Trmm
//   X := tril(L)^T,
//   X := tril(L)^H,
//   X := trilu(L)^T, or
//   X := trilu(L)^H
template<typename T>
inline void
TrmmLLT
( Orientation orientation, 
  UnitOrNonUnit diag,
  T alpha, 
  const DistMatrix<T>& L,
        DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmLLT");
#endif
    // TODO: Come up with a better routing mechanism
    if( L.Height() > 5*X.Width() )
        TrmmLLTA( orientation, diag, alpha, L, X );
    else
        TrmmLLTC( orientation, diag, alpha, L, X );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
