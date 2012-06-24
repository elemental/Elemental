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
HemmRUA
( T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::HemmRUA");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
#endif
    const Grid& g = A.Grid();

    DistMatrix<T>
        BT(g),  B0(g),
        BB(g),  B1(g),
                B2(g);
    DistMatrix<T>
        CT(g),  C0(g),
        CB(g),  C1(g),
                C2(g);

    DistMatrix<T,MR,  STAR> B1Adj_MR_STAR(g);
    DistMatrix<T,VC,  STAR> B1Adj_VC_STAR(g);
    DistMatrix<T,STAR,MC  > B1_STAR_MC(g);
    DistMatrix<T,MC,  STAR> Z1Adj_MC_STAR(g);
    DistMatrix<T,MR,  STAR> Z1Adj_MR_STAR(g);
    DistMatrix<T,MR,  MC  > Z1Adj_MR_MC(g);
    DistMatrix<T> Z1Adj(g);

    B1Adj_MR_STAR.AlignWith( A );
    B1Adj_VC_STAR.AlignWith( A );
    B1_STAR_MC.AlignWith( A );
    Z1Adj_MC_STAR.AlignWith( A );
    Z1Adj_MR_STAR.AlignWith( A );

    Matrix<T> Z1Local;

    Scale( beta, C );
    LockedPartitionDown
    ( B, BT,
         BB, 0 );
    PartitionDown
    ( C, CT,
         CB, 0 );
    while( CT.Height() < C.Height() )
    {
        LockedRepartitionDown
        ( BT,  B0, 
         /**/ /**/
               B1,
          BB,  B2 );

        RepartitionDown
        ( CT,  C0,
         /**/ /**/
               C1,
          CB,  C2 );

        Z1Adj_MR_MC.AlignWith( C1 );
        Zeros( C1.Width(), C1.Height(), Z1Adj_MC_STAR );
        Zeros( C1.Width(), C1.Height(), Z1Adj_MR_STAR );
        //--------------------------------------------------------------------//
        B1Adj_MR_STAR.AdjointFrom( B1 );
        B1Adj_VC_STAR = B1Adj_MR_STAR;
        B1_STAR_MC.AdjointFrom( B1Adj_VC_STAR );
        LocalSymmetricAccumulateRU
        ( ADJOINT, alpha, A, B1_STAR_MC, B1Adj_MR_STAR, 
          Z1Adj_MC_STAR, Z1Adj_MR_STAR );

        Z1Adj.SumScatterFrom( Z1Adj_MC_STAR );
        Z1Adj_MR_MC = Z1Adj;
        Z1Adj_MR_MC.SumScatterUpdate( (T)1, Z1Adj_MR_STAR );
        Adjoint( Z1Adj_MR_MC.LockedLocalMatrix(), Z1Local );
        Axpy( (T)1, Z1Local, C1.LocalMatrix() );
        //--------------------------------------------------------------------//
        Z1Adj_MR_MC.FreeAlignments();

        SlideLockedPartitionDown
        ( BT,  B0,
               B1,
         /**/ /**/
          BB,  B2 );

        SlidePartitionDown
        ( CT,  C0,
               C1,
         /**/ /**/
          CB,  C2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
HemmRUC
( T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::HemmRUC");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error("{A,B,C} must be distributed on the same grid");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  AColPan(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),  ARowPan(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<T> BL(g), BR(g),
                  B0(g), B1(g), B2(g);
    DistMatrix<T> CL(g), CR(g),
                  C0(g), C1(g), C2(g),
                  CLeft(g), CRight(g);

    // Temporary distributions
    DistMatrix<T,MC,STAR> B1_MC_STAR(g);
    DistMatrix<T,VR,  STAR> AColPan_VR_STAR(g);
    DistMatrix<T,STAR,MR  > AColPanAdj_STAR_MR(g);
    DistMatrix<T,MR,  STAR> ARowPanAdj_MR_STAR(g);

    B1_MC_STAR.AlignWith( C );

    // Start the algorithm
    Scale( beta, C );
    LockedPartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionRight( B, BL, BR, 0 );
    PartitionRight( C, CL, CR, 0 );
    while( CR.Width() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionRight
        ( BL, /**/ BR,
          B0, /**/ B1, B2 );

        RepartitionRight
        ( CL, /**/ CR,
          C0, /**/ C1, C2 );

        ARowPan.LockedView1x2( A11, A12 );
        AColPan.LockedView2x1
        ( A01,
          A11 );

        CLeft.View1x2( C0, C1 );
        CRight.View1x2( C1, C2 );

        AColPan_VR_STAR.AlignWith( CLeft );
        AColPanAdj_STAR_MR.AlignWith( CLeft );
        ARowPanAdj_MR_STAR.AlignWith( CRight );
        //--------------------------------------------------------------------//
        B1_MC_STAR = B1;

        AColPan_VR_STAR = AColPan;
        AColPanAdj_STAR_MR.AdjointFrom( AColPan_VR_STAR );
        ARowPanAdj_MR_STAR.AdjointFrom( ARowPan );
        MakeTrapezoidal( LEFT,  LOWER,  0, ARowPanAdj_MR_STAR );
        MakeTrapezoidal( RIGHT, LOWER, -1, AColPanAdj_STAR_MR );

        LocalGemm
        ( NORMAL, ADJOINT, 
          alpha, B1_MC_STAR, ARowPanAdj_MR_STAR, (T)1, CRight );

        LocalGemm
        ( NORMAL, NORMAL,
          alpha, B1_MC_STAR, AColPanAdj_STAR_MR, (T)1, CLeft );
        //--------------------------------------------------------------------//
        AColPan_VR_STAR.FreeAlignments();
        AColPanAdj_STAR_MR.FreeAlignments();
        ARowPanAdj_MR_STAR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );

        SlidePartitionRight
        ( CL,     /**/ CR,
          C0, C1, /**/ C2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalSymmetricAccumulateRU
( Orientation orientation, T alpha,
  const DistMatrix<T,MC,  MR  >& A,
  const DistMatrix<T,STAR,MC  >& B_STAR_MC,
  const DistMatrix<T,MR,  STAR>& BAdjOrTrans_MR_STAR,
        DistMatrix<T,MC,  STAR>& ZAdjOrTrans_MC_STAR,
        DistMatrix<T,MR,  STAR>& ZAdjOrTrans_MR_STAR )
{
#ifndef RELEASE
    PushCallStack("internal::LocalSymmetricAccumulateRU");
    if( A.Grid() != B_STAR_MC.Grid() ||
        B_STAR_MC.Grid() != BAdjOrTrans_MR_STAR.Grid() ||
        BAdjOrTrans_MR_STAR.Grid() != ZAdjOrTrans_MC_STAR.Grid() ||
        ZAdjOrTrans_MC_STAR.Grid() != ZAdjOrTrans_MR_STAR.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
    if( A.Height() != A.Width() ||
        A.Height() != B_STAR_MC.Width() ||
        A.Height() != BAdjOrTrans_MR_STAR.Height() ||
        A.Height() != ZAdjOrTrans_MC_STAR.Height() ||
        A.Height() != ZAdjOrTrans_MR_STAR.Height() ||
        B_STAR_MC.Height() != BAdjOrTrans_MR_STAR.Width() ||
        BAdjOrTrans_MR_STAR.Width() != ZAdjOrTrans_MC_STAR.Width() ||
        ZAdjOrTrans_MC_STAR.Width() != ZAdjOrTrans_MR_STAR.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalSymmetricAccumulateRU: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  B[* ,MC] ~ " << B_STAR_MC.Height() << " x "
                               << B_STAR_MC.Width() << "\n"
            << "  B^H/T[MR,* ] ~ " << BAdjOrTrans_MR_STAR.Height() << " x "
                                   << BAdjOrTrans_MR_STAR.Width() << "\n"
            << "  Z^H/T[MC,* ] ~ " << ZAdjOrTrans_MC_STAR.Height() << " x "
                                   << ZAdjOrTrans_MC_STAR.Width() << "\n"
            << "  Z^H/T[MR,* ] ~ " << ZAdjOrTrans_MR_STAR.Height() << " x "
                                   << ZAdjOrTrans_MR_STAR.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( B_STAR_MC.RowAlignment() != A.ColAlignment() ||
        BAdjOrTrans_MR_STAR.ColAlignment() != A.RowAlignment() ||
        ZAdjOrTrans_MC_STAR.ColAlignment() != A.ColAlignment() ||
        ZAdjOrTrans_MR_STAR.ColAlignment() != A.RowAlignment() )
        throw std::logic_error("Partial matrix distributions are misaligned");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<T> D11(g);

    DistMatrix<T,STAR,MC>
        BL_STAR_MC(g), BR_STAR_MC(g),
        B0_STAR_MC(g), B1_STAR_MC(g), B2_STAR_MC(g);

    DistMatrix<T,MR,STAR>
        BTAdjOrTrans_MR_STAR(g),  B0AdjOrTrans_MR_STAR(g),
        BBAdjOrTrans_MR_STAR(g),  B1AdjOrTrans_MR_STAR(g),
                                  B2AdjOrTrans_MR_STAR(g);

    DistMatrix<T,MC,STAR>
        ZTAdjOrTrans_MC_STAR(g),  Z0AdjOrTrans_MC_STAR(g),
        ZBAdjOrTrans_MC_STAR(g),  Z1AdjOrTrans_MC_STAR(g),
                                  Z2AdjOrTrans_MC_STAR(g);

    DistMatrix<T,MR,STAR>
        ZBAdjOrTrans_MR_STAR(g),  Z0AdjOrTrans_MR_STAR(g),
        ZTAdjOrTrans_MR_STAR(g),  Z1AdjOrTrans_MR_STAR(g),
                                  Z2AdjOrTrans_MR_STAR(g);

    const int ratio = std::max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*Blocksize() );

    LockedPartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionRight( B_STAR_MC,  BL_STAR_MC, BR_STAR_MC, 0 );
    LockedPartitionDown
    ( BAdjOrTrans_MR_STAR, BTAdjOrTrans_MR_STAR,
                           BBAdjOrTrans_MR_STAR, 0 );
    PartitionDown
    ( ZAdjOrTrans_MC_STAR, ZTAdjOrTrans_MC_STAR,
                           ZBAdjOrTrans_MC_STAR, 0 );
    PartitionDown
    ( ZAdjOrTrans_MR_STAR, ZTAdjOrTrans_MR_STAR,
                           ZBAdjOrTrans_MR_STAR, 0 );
    while( ATL.Height() < A.Height() )
    {
        LockedRepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionRight
        ( BL_STAR_MC, /**/ BR_STAR_MC,
          B0_STAR_MC, /**/ B1_STAR_MC, B2_STAR_MC );

        LockedRepartitionDown
        ( BTAdjOrTrans_MR_STAR,  B0AdjOrTrans_MR_STAR,
         /********************/ /********************/
                                 B1AdjOrTrans_MR_STAR,
          BBAdjOrTrans_MR_STAR,  B2AdjOrTrans_MR_STAR );

        RepartitionDown
        ( ZTAdjOrTrans_MC_STAR,  Z0AdjOrTrans_MC_STAR,
         /********************/ /********************/
                                 Z1AdjOrTrans_MC_STAR,
          ZBAdjOrTrans_MC_STAR,  Z2AdjOrTrans_MC_STAR );

        RepartitionDown
        ( ZTAdjOrTrans_MR_STAR,  Z0AdjOrTrans_MR_STAR,
         /********************/ /********************/
                                 Z1AdjOrTrans_MR_STAR,
          ZBAdjOrTrans_MR_STAR,  Z2AdjOrTrans_MR_STAR );

        D11.AlignWith( A11 );
        //--------------------------------------------------------------------//
        D11 = A11;
        MakeTrapezoidal( LEFT, UPPER, 0, D11 );
        LocalGemm
        ( orientation, orientation,
          alpha, D11, B1_STAR_MC, (T)1, Z1AdjOrTrans_MR_STAR );
        MakeTrapezoidal( LEFT, UPPER, 1, D11 );

        LocalGemm
        ( NORMAL, NORMAL, alpha, D11, B1AdjOrTrans_MR_STAR, 
          (T)1, Z1AdjOrTrans_MC_STAR );

        LocalGemm
        ( orientation, orientation,
          alpha, A12, B1_STAR_MC, (T)1, Z2AdjOrTrans_MR_STAR );

        LocalGemm
        ( NORMAL, NORMAL, alpha, A12, B2AdjOrTrans_MR_STAR, 
          (T)1, Z1AdjOrTrans_MC_STAR );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionRight
        ( BL_STAR_MC,             /**/ BR_STAR_MC,
          B0_STAR_MC, B1_STAR_MC, /**/ B2_STAR_MC );

        SlideLockedPartitionDown
        ( BTAdjOrTrans_MR_STAR,  B0AdjOrTrans_MR_STAR,
                                 B1AdjOrTrans_MR_STAR,
         /********************/ /********************/
          BBAdjOrTrans_MR_STAR,  B2AdjOrTrans_MR_STAR );

        SlidePartitionDown
        ( ZTAdjOrTrans_MC_STAR,  Z0AdjOrTrans_MC_STAR,
                                 Z1AdjOrTrans_MC_STAR,
         /********************/ /********************/
          ZBAdjOrTrans_MC_STAR,  Z2AdjOrTrans_MC_STAR );       
        
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

template<typename T>
inline void
HemmRU
( T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::HemmRU");
#endif
    // TODO: Come up with a better routing mechanism
    if( A.Height() > 5*B.Height() )
        HemmRUA( alpha, A, B, beta, C );
    else
        HemmRUC( alpha, A, B, beta, C );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
