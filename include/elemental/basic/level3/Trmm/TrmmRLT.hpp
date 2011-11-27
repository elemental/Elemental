/*
   Copyright (c) 2009-2011, Jack Poulson
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

// Right Lower Adjoint/Transpose (Non)Unit Trmm
//   X := X tril(L)^T,
//   X := X tril(L)^H,
//   X := X trilu(L)^T, or
//   X := X trilu(L)^H
template<typename T>
inline void
elemental::basic::internal::TrmmRLT
( Orientation orientation, 
  Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("basic::internal::TrmmRLT");
#endif
    // TODO: Come up with a better routing mechanism
    if( L.Height() > 5*X.Height() )
        basic::internal::TrmmRLTA( orientation, diagonal, alpha, L, X );
    else
        basic::internal::TrmmRLTC( orientation, diagonal, alpha, L, X );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::internal::TrmmRLTA
( Orientation orientation, Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("basic::internal::TrmmRLTA");
    if( L.Grid() != X.Grid() )
        throw std::logic_error("{L,X} must be distributed over the same grid");
#endif
    const Grid& g = L.Grid();

    DistMatrix<T,MC,MR>
        XT(g),  X0(g),
        XB(g),  X1(g),
                X2(g);

    DistMatrix<T,MR,  STAR> X1AdjOrTrans_MR_STAR(g);
    DistMatrix<T,MC,  STAR> Z1AdjOrTrans_MC_STAR(g);
    DistMatrix<T,MC,  MR  > Z1AdjOrTrans(g);
    DistMatrix<T,MR,  MC  > Z1AdjOrTrans_MR_MC(g);

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

        X1AdjOrTrans_MR_STAR.AlignWith( L );
        Z1AdjOrTrans_MC_STAR.AlignWith( L );
        Z1AdjOrTrans_MR_MC.AlignWith( X1 );
        Z1AdjOrTrans_MC_STAR.ResizeTo( X1.Width(), X1.Height() );
        //--------------------------------------------------------------------//
        if( orientation == ADJOINT )
            X1AdjOrTrans_MR_STAR.AdjointFrom( X1 );
        else
            X1AdjOrTrans_MR_STAR.TransposeFrom( X1 );
        Z1AdjOrTrans_MC_STAR.SetToZero();
        basic::internal::LocalTrmmAccumulateRLT
        ( orientation, diagonal, alpha, 
          L, X1AdjOrTrans_MR_STAR, Z1AdjOrTrans_MC_STAR );

        Z1AdjOrTrans.SumScatterFrom( Z1AdjOrTrans_MC_STAR );
        Z1AdjOrTrans_MR_MC = Z1AdjOrTrans;
        if( orientation == ADJOINT )
            basic::Adjoint
            ( Z1AdjOrTrans_MR_MC.LocalMatrix(), X1.LocalMatrix() );
        else
            basic::Transpose
            ( Z1AdjOrTrans_MR_MC.LocalMatrix(), X1.LocalMatrix() );
        //--------------------------------------------------------------------//
        X1AdjOrTrans_MR_STAR.FreeAlignments();
        Z1AdjOrTrans_MC_STAR.FreeAlignments();
        Z1AdjOrTrans_MR_MC.FreeAlignments();

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
elemental::basic::internal::TrmmRLTC
( Orientation orientation, 
  Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("basic::internal::TrmmRLTC");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( orientation == NORMAL )
        throw std::logic_error("TrmmRLTC expects an Adjoint/Transpose option");
    if( L.Height() != L.Width() || X.Width() != L.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal TrmmRLTC: \n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<T,MC,MR> XL(g), XR(g),
                        X0(g), X1(g), X2(g);

    // Temporary distributions
    DistMatrix<T,STAR,MR  > L10_STAR_MR(g);
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<T,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<T,MC,  STAR> D1_MC_STAR(g);

    // Start the algorithm
    basic::Scal( alpha, X );
    LockedPartitionUpDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionLeft( X, XL, XR, 0 );
    while( XL.Width() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        RepartitionLeft
        ( XL,     /**/ XR,
          X0, X1, /**/ X2 );

        L10_STAR_MR.AlignWith( X0 );
        D1_MC_STAR.AlignWith( X1 );
        D1_MC_STAR.ResizeTo( X1.Height(), X1.Width() );
        //--------------------------------------------------------------------//
        X1_VC_STAR = X1;
        L11_STAR_STAR = L11;
        basic::internal::LocalTrmm
        ( RIGHT, LOWER, orientation, diagonal, 
          (T)1, L11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
 
        L10_STAR_MR = L10;
        basic::internal::LocalGemm
        ( NORMAL, orientation, (T)1, X0, L10_STAR_MR, (T)0, D1_MC_STAR );
        X1.SumScatterUpdate( (T)1, D1_MC_STAR );
       //--------------------------------------------------------------------//
        L10_STAR_MR.FreeAlignments();
        D1_MC_STAR.FreeAlignments();

        SlideLockedPartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        SlidePartitionLeft
        ( XL, /**/     XR,
          X0, /**/ X1, X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::internal::LocalTrmmAccumulateRLT
( Orientation orientation, Diagonal diagonal, T alpha,
  const DistMatrix<T,MC,MR  >& L,
  const DistMatrix<T,MR,STAR>& XAdjOrTrans_MR_STAR,
        DistMatrix<T,MC,STAR>& ZAdjOrTrans_MC_STAR )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalTrmmAccumulateRLT");
    if( L.Grid() != XAdjOrTrans_MR_STAR.Grid() ||
        XAdjOrTrans_MR_STAR.Grid() != ZAdjOrTrans_MC_STAR.Grid() )
        throw std::logic_error
        ("{L,X,Z} must be distributed over the same grid");
    if( L.Height() != L.Width() ||
        L.Height() != XAdjOrTrans_MR_STAR.Height() ||
        L.Height() != ZAdjOrTrans_MC_STAR.Height() ||
        XAdjOrTrans_MR_STAR.Width() != ZAdjOrTrans_MC_STAR.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrmmAccumulateRLT: \n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X^H/T[MR,* ] ~ " << XAdjOrTrans_MR_STAR.Height() << " x "
                                   << XAdjOrTrans_MR_STAR.Width() << "\n"
            << "  Z^H/T[MC,* ] ~ " << ZAdjOrTrans_MC_STAR.Height() << " x "
                                   << ZAdjOrTrans_MC_STAR.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( XAdjOrTrans_MR_STAR.ColAlignment() != L.RowAlignment() ||
        ZAdjOrTrans_MC_STAR.ColAlignment() != L.ColAlignment() )
        throw std::logic_error("Partial matrix distributions are misaligned");
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<T,MC,MR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<T,MC,MR> D11(g);

    DistMatrix<T,MR,STAR>
        XTAdjOrTrans_MR_STAR(g),  X0AdjOrTrans_MR_STAR(g),
        XBAdjOrTrans_MR_STAR(g),  X1AdjOrTrans_MR_STAR(g),
                                   X2AdjOrTrans_MR_STAR(g);

    DistMatrix<T,MC,STAR>
        ZTAdjOrTrans_MC_STAR(g),  Z0AdjOrTrans_MC_STAR(g),
        ZBAdjOrTrans_MC_STAR(g),  Z1AdjOrTrans_MC_STAR(g),
                                   Z2AdjOrTrans_MC_STAR(g);

    const int ratio = std::max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );

    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    LockedPartitionDown
    ( XAdjOrTrans_MR_STAR, XTAdjOrTrans_MR_STAR,
                           XBAdjOrTrans_MR_STAR, 0 );
    PartitionDown
    ( ZAdjOrTrans_MC_STAR, ZTAdjOrTrans_MC_STAR,
                           ZBAdjOrTrans_MC_STAR, 0 );
    while( LTL.Height() < L.Height() )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        LockedRepartitionDown
        ( XTAdjOrTrans_MR_STAR,  X0AdjOrTrans_MR_STAR,
         /********************/ /********************/
                                 X1AdjOrTrans_MR_STAR,
          XBAdjOrTrans_MR_STAR,  X2AdjOrTrans_MR_STAR );

        RepartitionDown
        ( ZTAdjOrTrans_MC_STAR,  Z0AdjOrTrans_MC_STAR,
         /********************/ /********************/
                                 Z1AdjOrTrans_MC_STAR,
          ZBAdjOrTrans_MC_STAR,  Z2AdjOrTrans_MC_STAR );

        D11.AlignWith( L11 );
        //--------------------------------------------------------------------//
        D11 = L11;
        D11.MakeTrapezoidal( LEFT, LOWER );
        if( diagonal == UNIT )
            SetDiagonalToOne( D11 );

        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, D11, X1AdjOrTrans_MR_STAR, 
          (T)1, Z1AdjOrTrans_MC_STAR );

        basic::internal::LocalGemm
        ( NORMAL, NORMAL, alpha, L21, X1AdjOrTrans_MR_STAR, 
          (T)1, Z2AdjOrTrans_MC_STAR );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlideLockedPartitionDown
        ( XTAdjOrTrans_MR_STAR,   X0AdjOrTrans_MR_STAR,
                                  X1AdjOrTrans_MR_STAR,
         /*********************/ /********************/
          XBAdjOrTrans_MR_STAR,   X2AdjOrTrans_MR_STAR );

        SlidePartitionDown
        ( ZTAdjOrTrans_MC_STAR,  Z0AdjOrTrans_MC_STAR,
                                 Z1AdjOrTrans_MC_STAR,
         /********************/ /********************/
          ZBAdjOrTrans_MC_STAR,  Z2AdjOrTrans_MC_STAR );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}
