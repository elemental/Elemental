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
TrmmRUTA
( Orientation orientation, UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& U,
                 DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmRUTA");
    if( U.Grid() != X.Grid() )
        throw std::logic_error("{U,X} must be distributed over the same grid");
#endif
    const Grid& g = U.Grid();

    DistMatrix<T>
        XT(g),  X0(g),
        XB(g),  X1(g),
                X2(g);

    DistMatrix<T,MR,  STAR> X1AdjOrTrans_MR_STAR(g);
    DistMatrix<T,MC,  STAR> Z1AdjOrTrans_MC_STAR(g);
    DistMatrix<T,MC,  MR  > Z1AdjOrTrans(g);
    DistMatrix<T,MR,  MC  > Z1AdjOrTrans_MR_MC(g);

    X1AdjOrTrans_MR_STAR.AlignWith( U );
    Z1AdjOrTrans_MC_STAR.AlignWith( U );

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

        Z1AdjOrTrans_MR_MC.AlignWith( X1 );
        Zeros( X1.Width(), X1.Height(), Z1AdjOrTrans_MC_STAR );
        //--------------------------------------------------------------------//
        if( orientation == ADJOINT )
            X1AdjOrTrans_MR_STAR.AdjointFrom( X1 );
        else
            X1AdjOrTrans_MR_STAR.TransposeFrom( X1 );
        LocalTrmmAccumulateRUT
        ( orientation, diag, alpha,
          U, X1AdjOrTrans_MR_STAR, Z1AdjOrTrans_MC_STAR );

        Z1AdjOrTrans.SumScatterFrom( Z1AdjOrTrans_MC_STAR );
        Z1AdjOrTrans_MR_MC = Z1AdjOrTrans;
        if( orientation == ADJOINT )
            Adjoint( Z1AdjOrTrans_MR_MC.LocalMatrix(), X1.LocalMatrix() );
        else
            Transpose( Z1AdjOrTrans_MR_MC.LocalMatrix(), X1.LocalMatrix() );
        //--------------------------------------------------------------------//
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
TrmmRUTC
( Orientation orientation, 
  UnitOrNonUnit diag,
  T alpha, const DistMatrix<T>& U,
                 DistMatrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("internal::TrmmRUTC");
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

    // Matrix views
    DistMatrix<T> 
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    DistMatrix<T> XL(g), XR(g),
                  X0(g), X1(g), X2(g);

    // Temporary distributions
    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<T,MR,  STAR> U12AdjOrTrans_MR_STAR(g);
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

        U12AdjOrTrans_MR_STAR.AlignWith( X2 );
        D1_MC_STAR.AlignWith( X1 );
        Zeros( X1.Height(), X1.Width(), D1_MC_STAR );
        //--------------------------------------------------------------------//
        X1_VC_STAR = X1;
        U11_STAR_STAR = U11;
        LocalTrmm
        ( RIGHT, UPPER, orientation, diag, T(1), U11_STAR_STAR, X1_VC_STAR );
        X1 = X1_VC_STAR;
 
        if( orientation == ADJOINT )
            U12AdjOrTrans_MR_STAR.AdjointFrom( U12 );
        else
            U12AdjOrTrans_MR_STAR.TransposeFrom( U12 );
        LocalGemm
        ( NORMAL, NORMAL, T(1), X2, U12AdjOrTrans_MR_STAR, T(0), D1_MC_STAR );
        X1.SumScatterUpdate( T(1), D1_MC_STAR );
        //--------------------------------------------------------------------//
        U12AdjOrTrans_MR_STAR.FreeAlignments();
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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalTrmmAccumulateRUT
( Orientation orientation, UnitOrNonUnit diag, T alpha,
  const DistMatrix<T>& U,
  const DistMatrix<T,MR,STAR>& XAdjOrTrans_MR_STAR,
        DistMatrix<T,MC,STAR>& ZAdjOrTrans_MC_STAR )
{
#ifndef RELEASE
    PushCallStack("internal::LocalTrmmAccumulateRUT");
    if( U.Grid() != XAdjOrTrans_MR_STAR.Grid() ||
        XAdjOrTrans_MR_STAR.Grid() != ZAdjOrTrans_MC_STAR.Grid() )
        throw std::logic_error
        ("{U,X,Z} must be distributed over the same grid");
    if( U.Height() != U.Width() ||
        U.Height() != XAdjOrTrans_MR_STAR.Height() ||
        U.Height() != ZAdjOrTrans_MC_STAR.Height() ||
        XAdjOrTrans_MR_STAR.Width() != ZAdjOrTrans_MC_STAR.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalTrmmAccumulateRUT: \n"
            << "  U ~ " << U.Height() << " x " << U.Width() << "\n"
            << "  X^H/T[MR,* ] ~ " << XAdjOrTrans_MR_STAR.Height() << " x "
                                   << XAdjOrTrans_MR_STAR.Width() << "\n"
            << "  Z^H/T[MC,* ] ~ " << ZAdjOrTrans_MC_STAR.Height() << " x "
                                   << ZAdjOrTrans_MC_STAR.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( XAdjOrTrans_MR_STAR.ColAlignment() != U.RowAlignment() ||
        ZAdjOrTrans_MC_STAR.ColAlignment() != U.ColAlignment() )
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
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    LockedPartitionDown
    ( XAdjOrTrans_MR_STAR, XTAdjOrTrans_MR_STAR,
                           XBAdjOrTrans_MR_STAR, 0 );
    PartitionDown
    ( ZAdjOrTrans_MC_STAR, ZTAdjOrTrans_MC_STAR,
                           ZBAdjOrTrans_MC_STAR, 0 );
    while( UTL.Height() < U.Height() )
    {
        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

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

        D11.AlignWith( U11 );
        //--------------------------------------------------------------------//
        D11 = U11;
        MakeTrapezoidal( LEFT, UPPER, 0, D11 );
        if( diag == UNIT )
            SetDiagonalToOne( D11 );

        LocalGemm
        ( NORMAL, NORMAL, alpha, D11, X1AdjOrTrans_MR_STAR,
          T(1), Z1AdjOrTrans_MC_STAR );

        LocalGemm
        ( NORMAL, NORMAL, alpha, U01, X1AdjOrTrans_MR_STAR,
          T(1), Z0AdjOrTrans_MC_STAR );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

        SlideLockedPartitionDown
        ( XTAdjOrTrans_MR_STAR,  X0AdjOrTrans_MR_STAR,
                                 X1AdjOrTrans_MR_STAR,
         /********************/ /********************/
          XBAdjOrTrans_MR_STAR,  X2AdjOrTrans_MR_STAR );

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
    PushCallStack("internal::TrmmRUT");
#endif
    // TODO: Come up with a better routing mechanism
    if( U.Height() > 5*X.Height() )
        TrmmRUTA( orientation, diag, alpha, U, X );
    else
        TrmmRUTC( orientation, diag, alpha, U, X );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
