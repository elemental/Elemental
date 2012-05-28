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

template<typename R>
inline void 
PanelBidiagU
( DistMatrix<R,MC,  MR  >& A, 
  DistMatrix<R,MC,  MR  >& X, 
  DistMatrix<R,MC,  MR  >& Y,
  DistMatrix<R,MC,  STAR>& AColPan_MC_STAR,
  DistMatrix<R,STAR,MR  >& ARowPan_STAR_MR )
{
    const int panelSize = X.Width();
#ifndef RELEASE
    PushCallStack("internal::PanelBidiagU");
    if( A.Grid() != X.Grid() || X.Grid() != Y.Grid() ||
        Y.Grid() != AColPan_MC_STAR.Grid() || 
        Y.Grid() != ARowPan_STAR_MR.Grid() )
        throw std::logic_error("Grids must match");
    if( A.Height() < A.Width() )
        throw std::logic_error("A must be at least as tall as it is wide");
    if( A.Height() != X.Height() )
        throw std::logic_error("A and X must be the same height");
    if( A.Width() != Y.Height() )
        throw std::logic_error("Y must be the same height as A's width");
    if( X.Height() < panelSize )
        throw std::logic_error("X must be a column panel");
    if( Y.Width() != panelSize )
        throw std::logic_error("Y is the wrong width");
    if( A.ColAlignment() != X.ColAlignment() || 
        A.RowAlignment() != X.RowAlignment() )
        throw std::logic_error("A and X must be aligned");
    if( A.ColAlignment() != Y.ColAlignment() ||
        A.RowAlignment() != Y.RowAlignment() )
        throw std::logic_error("A and Y must be aligned");
#endif
    const Grid& g = A.Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int p = g.Size();

    // Matrix views 
    DistMatrix<R,MC,MR> 
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  aB1(g), AB2(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),  alpha12L(g), a12R(g),
                         A20(g), a21(g),     A22(g),  A2L(g);
    DistMatrix<R,MC,MR>
        XTL(g), XTR(g),  X00(g), x01(g),   X02(g), 
        XBL(g), XBR(g),  x10(g), chi11(g), x12(g), 
                         X20(g), x21(g),   X22(g);
    DistMatrix<R,MC,MR>
        YTL(g), YTR(g),  Y00(g), y01(g),   Y02(g),
        YBL(g), YBR(g),  y10(g), psi11(g), y12(g),
                         Y20(g), y21(g),   Y22(g),  Y2L(g);

    DistMatrix<R,MD,STAR> d(g), dT(g), d0(g), 
                                dB(g), delta1(g),
                                       d2(g);
    DistMatrix<R,MD,STAR> e(g), eT(g), e0(g),
                                eB(g), epsilon1(g),
                                       e2(g);
    DistMatrix<R,MC,STAR> aB1_MC_STAR(g);
    DistMatrix<R,STAR,MR> a12_STAR_MR(g);

    // Temporary distributions
    DistMatrix<R,MR,  STAR> a01_MR_STAR(g);
    DistMatrix<R,STAR,MR  > a10_STAR_MR(g);
    DistMatrix<R,STAR,MC  > a12_STAR_MC(g);
    DistMatrix<R,STAR,MC  > x10_STAR_MC(g);
    DistMatrix<R,STAR,MR  > y10_STAR_MR(g);
    DistMatrix<R,MC,  STAR> uB1_MC_STAR(g);
    DistMatrix<R,MR,  MC  > z01_MR_MC(g);
    DistMatrix<R,MC,  STAR> z01_MC_STAR(g);
    DistMatrix<R,MR,  STAR> z01_MR_STAR(g);
    DistMatrix<R,MR,  MC  > z21_MR_MC(g);
    DistMatrix<R,MC,  STAR> z21_MC_STAR(g);
    DistMatrix<R,MR,  STAR> z21_MR_STAR(g);
    DistMatrix<R,MC,  MR  > q21(g);
    DistMatrix<R,MR,  MC  > q21_MR_MC(g);
    DistMatrix<R,MC,  STAR> q21_MC_STAR(g);
    DistMatrix<R,MR,  STAR> q21_MR_STAR(g);
    DistMatrix<R,MC,  MR  > s01(g);
    DistMatrix<R,MC,  STAR> s01_MC_STAR(g);
    DistMatrix<R,MR,  STAR> s01_MR_STAR(g);
    DistMatrix<R,MC,  STAR> s21_MC_STAR(g);
    DistMatrix<R,MR,  STAR> sB1_MR_STAR(g);

    d.AlignWithDiagonal( A, 0 );
    e.AlignWithDiagonal( A, 1 );
    d.ResizeTo( panelSize, 1 );
    e.ResizeTo( panelSize, 1 );

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDownDiagonal
    ( X, XTL, XTR,
         XBL, XBR, 0 );
    PartitionDownDiagonal
    ( Y, YTL, YTR,
         YBL, YBR, 0 );
    PartitionDown
    ( d, dT,
         dB, 0 );
    PartitionDown
    ( e, eT,
         eB, 0 );
    PushBlocksizeStack( 1 );
    while( ATL.Width() < panelSize )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        RepartitionDownDiagonal
        ( XTL, /**/ XTR,  X00, /**/ x01,   X02,
         /*************/ /********************/
               /**/       x10, /**/ chi11, x12,
          XBL, /**/ XBR,  X20, /**/ x21,   X22 );
        
        RepartitionDownDiagonal
        ( YTL, /**/ YTR,  Y00, /**/ y01,   Y02,
         /*************/ /********************/
               /**/       y10, /**/ psi11, y12,
          YBL, /**/ YBR,  Y20, /**/ y21,   Y22 );

        RepartitionDown
        ( dT,  d0,
         /**/ /******/
               delta1,
          dB,  d2 );

        RepartitionDown
        ( eT,  e0,
         /**/ /********/
               epsilon1,
          eB,  e2 );

        PartitionRight( ABR, aB1, AB2, 1 );
        PartitionRight( a12, alpha12L, a12R, 1 );

        A2L.View1x2( A20, a21 );
        Y2L.View1x2( Y20, y21 );

        a12_STAR_MR.View
        ( ARowPan_STAR_MR, ATL.Height(), ATL.Width()+1, 1, a12.Width() );
        aB1_MC_STAR.View
        ( AColPan_MC_STAR, ATL.Height(), ATL.Width(), ABR.Height(), 1 );

        // Main alignments
        a01_MR_STAR.AlignWith( ABL );
        a10_STAR_MR.AlignWith( Y20 );
        a12_STAR_MC.AlignWith( Y2L );
        x10_STAR_MC.AlignWith( A02 );
        y10_STAR_MR.AlignWith( ABL );

        // Auxilliary alignments
        uB1_MC_STAR.AlignWith( ABL );
        z01_MC_STAR.AlignWith( A02 );
        z01_MR_STAR.AlignWith( ABL );
        z21_MC_STAR.AlignWith( Y20 );
        z21_MR_STAR.AlignWith( AB2 );
        q21.AlignWith( y21 );
        q21_MR_MC.AlignWith( a12 );
        q21_MC_STAR.AlignWith( Y20 );
        q21_MR_STAR.AlignWith( A02 );
        s01_MC_STAR.AlignWith( A02 );
        s01_MR_STAR.AlignWith( X20 );
        s21_MC_STAR.AlignWith( A22 );
        sB1_MR_STAR.AlignWith( Y2L );

        // Auxilliary resizes
        uB1_MC_STAR.ResizeTo( ABL.Height(), 1 );
        z01_MR_STAR.ResizeTo( A00.Width(), 1 );
        z21_MC_STAR.ResizeTo( A22.Width(), 1 );
        z21_MR_STAR.ResizeTo( A22.Width(), 1 );
        q21_MC_STAR.ResizeTo( A22.Width(), 1 );
        q21_MR_STAR.ResizeTo( A22.Width(), 1 );
        s01_MC_STAR.ResizeTo( A00.Height(), 1 );
        s21_MC_STAR.ResizeTo( A22.Height(), 1 );
        sB1_MR_STAR.ResizeTo( Y2L.Width(), 1 );

        const bool thisIsMyRow = ( g.Row() == alpha11.ColAlignment() );
        const bool thisIsMyCol = ( g.Col() == alpha11.RowAlignment() );
        const bool nextIsMyCol = ( g.Col() == a12.RowAlignment() ) ;
        const bool firstIteration = ( ATL.Height() == 0 );
        //--------------------------------------------------------------------//

        // Update the current column of A:
        //   aB1 := aB1 - ABL y10^T - XBL a01
        if( !firstIteration )
        {
            y10_STAR_MR = y10;
            a01_MR_STAR = a01;
            // uB1[MC,* ] := ABL[MC,MR] y10^T[MR,* ]
            Gemv
            ( NORMAL, 
              (R)1, ABL.LocalMatrix(), y10_STAR_MR.LocalMatrix(), 
              (R)0, uB1_MC_STAR.LocalMatrix() );
            // uB1[MC,* ] := uB1[MC,* ] + XBL[MC,MR] a01[MR,* ]
            //             = ABL[MC,MR] y10^T[MR,* ] + XBL[MC,MR] a01[MR,* ]
            Gemv
            ( NORMAL,
              (R)1, XBL.LocalMatrix(), a01_MR_STAR.LocalMatrix(),
              (R)1, uB1_MC_STAR.LocalMatrix() );
            // Sum the partial contributions and subtract from aB1
            aB1.SumScatterUpdate( (R)-1, uB1_MC_STAR );
        }

        // Find tauQ, u, and delta such that
        //     I - tauQ | 1 | | 1, u^T | | alpha11 | = | delta |
        //              | u |            |   a21   |   |   0   |
        R tauQ = 0;
        if( thisIsMyCol )
        {
            tauQ = ColReflector( alpha11, a21 );
            if( thisIsMyRow )
            {
                // Store delta and force | alpha11 | = | 1 |
                //                       |   a21   |   | u |
                delta1.SetLocalEntry(0,0,alpha11.GetLocalEntry(0,0));
                alpha11.SetLocalEntry(0,0,(R)1);
            }
        }

        //
        // y21 := tauQ ( AB2^T aB1 - A02^T XBL^T aB1 - Y20 ABL^T aB1 )
        //
        aB1_MC_STAR = aB1;
        // z01[MR,* ] := ABL^T[MR,MC] aB1[MC,* ]
        Gemv
        ( TRANSPOSE, 
          (R)1, ABL.LocalMatrix(), aB1_MC_STAR.LocalMatrix(), 
          (R)0, z01_MR_STAR.LocalMatrix() );
        // z21[MR,* ] := AB2^T[MR,MC] aB1[MC,* ]
        Gemv
        ( TRANSPOSE, 
          (R)1, AB2.LocalMatrix(), aB1_MC_STAR.LocalMatrix(), 
          (R)0, z21_MR_STAR.LocalMatrix() );
        // Sum the partial contributions to z01[MR,* ]
        z01_MR_STAR.SumOverCol();
        // z21[MC,* ] := Y20[MC,MR] z01[MR,* ] = Y20[MC,MR] (ABL^T aB1)[MR,* ]
        Gemv 
        ( NORMAL, 
          (R)1, Y20.LocalMatrix(), z01_MR_STAR.LocalMatrix(),
          (R)0, z21_MC_STAR.LocalMatrix() );
        // z01[MR,* ] := XBL^T[MR,MC] aB1[MC,* ]
        Gemv
        ( TRANSPOSE,
          (R)1, XBL.LocalMatrix(), aB1_MC_STAR.LocalMatrix(),
          (R)0, z01_MR_STAR.LocalMatrix() );
        // Sum the partial contributions to z01[MR,* ] and scatter the result
        z01_MR_MC.SumScatterFrom( z01_MR_STAR );
        // Redistribute the scattered summation 
        z01_MC_STAR = z01_MR_MC;
        // z21[MR,* ] := z21[MR,* ] - A02^T[MR,MC] z01[MC,* ] 
        //             = AB2^T[MR,MC] aB1[MC,* ] - 
        //               A02^T[MR,MC] (XBL^T aB1)[MC,* ]
        Gemv
        ( TRANSPOSE,
          (R)-1, A02.LocalMatrix(), z01_MC_STAR.LocalMatrix(),
          (R)1, z21_MR_STAR.LocalMatrix() );
        // Sum the partial contributions to z21[MR,* ] and scatter the result
        z21_MR_MC.SumScatterFrom( z21_MR_STAR );
        // Redistribute (and rename) the scattered summation
        y21 = z21_MR_MC;
        // Substract z21 = Y20 ABL^T aB1 from y21
        y21.SumScatterUpdate( (R)-1, z21_MC_STAR );
        if( thisIsMyCol )
            Scal( tauQ, y21 );

        // 
        // y21 := y21 + Y20 a10^T
        //
        a10_STAR_MR = a10;
        x10_STAR_MC = x10;
        // q21[MC,* ] := Y20[MC,MR] a10^T[MR,* ]
        Gemv
        ( NORMAL, 
          (R)1, Y20.LocalMatrix(), a10_STAR_MR.LocalMatrix(),
          (R)0, q21_MC_STAR.LocalMatrix() );
        // Sum the partial contributions
        q21.SumScatterFrom( q21_MC_STAR );
        if( thisIsMyCol )
            Axpy( (R)1, y21, q21 );

        //
        // a12 := a12 - a10 Y20^T - x10 A02
        //
        q21_MR_MC = q21;
        // q21[MR,* ] := A02^T[MR,MC] x10^T[MC,* ]
        Gemv
        ( TRANSPOSE,
          (R)1, A02.LocalMatrix(), x10_STAR_MC.LocalMatrix(),
          (R)0, q21_MR_STAR.LocalMatrix() );
        // Sum the partial contributions onto q21[MR,MC] = (Y20 a10^T)[MR,MC]
        q21_MR_MC.SumScatterUpdate( (R)1, q21_MR_STAR );
        // a12 := a12 - q21^T
        if( thisIsMyRow )
        {
            const int localWidth = a12.LocalWidth();
            R* a12Buffer = a12.LocalBuffer();
            const R* q21Buffer = q21_MR_MC.LockedLocalBuffer();
            const int a12LDim = a12.LocalLDim();
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                a12Buffer[jLocal*a12LDim] -= q21Buffer[jLocal];
        }

        // Find tauP, v, and epsilon such that
        //     I - tauP | 1 | | 1, v^T | | alpha12L | = | epsilon |
        //              | v |            |  a12R^T  |   |    0    |
        R tauP = 0;
        if( thisIsMyRow )
        {
            tauP = RowReflector( alpha12L, a12R );
            if( nextIsMyCol )
            {
                // Store epsilon and force | alpha12L | = | 1 |
                //                         |  a21R^T  |   | v |
                epsilon1.SetLocalEntry(0,0,alpha12L.GetLocalEntry(0,0));
                alpha12L.SetLocalEntry(0,0,(R)1);
            }
        }
        mpi::Broadcast( &tauP, 1, alpha11.ColAlignment(), g.ColComm() );

        //
        // x21 := tauP ( A22 a12^T - A2L Y2L^T a12^T - X20 A02 a12^T )
        //
        a12_STAR_MR = a12;
        a12_STAR_MC = a12;
        // s21[MC,* ] := A22[MC,MR] a12^T[MR,* ]
        Gemv
        ( NORMAL,
          (R)1, A22.LocalMatrix(), a12_STAR_MR.LocalMatrix(),
          (R)0, s21_MC_STAR.LocalMatrix() );
        // sB1[MR,* ] := Y2L^T[MR,MC] a12^T[MC,* ]
        Gemv
        ( TRANSPOSE,
          (R)1, Y2L.LocalMatrix(), a12_STAR_MC.LocalMatrix(),
          (R)0, sB1_MR_STAR.LocalMatrix() );
        // Sum the partial contributions
        sB1_MR_STAR.SumOverCol(); 
        // s21[MC,* ] := s21[MC,* ] - A2L[MC,MR] sB1[MR,* ]
        //             = A22[MC,MR] a12^T[MR,* ] - A2L[MC,MR] sB1[MR,* ]
        // (still needs to be summed within each process row)
        Gemv
        ( NORMAL, 
          (R)-1, A2L.LocalMatrix(), sB1_MR_STAR.LocalMatrix(),
          (R)1,  s21_MC_STAR.LocalMatrix() );
        // s01[MC,* ] := A02[MC,MR] a12^T[MR,* ]
        Gemv
        ( NORMAL,
          (R)1, A02.LocalMatrix(), a12_STAR_MR.LocalMatrix(),
          (R)0, s01_MC_STAR.LocalMatrix() );
        // Sum the partial contributions and then redistribute
        s01.SumScatterFrom( s01_MC_STAR ); // TODO: SumScatter to [VC,* ]?
        s01_MR_STAR = s01;
        // s21[MC,* ] := s21[MC,* ] - X20[MC,MR] s01[MR,* ]
        //             = A22[MC,MR] a12^T[MR,* ] - A2L[MC,MR] sB1[MR,* ]
        //                                       - X20[MC,MR] s01[MR,* ]
        Gemv
        ( NORMAL,
          (R)-1, X20.LocalMatrix(), s01_MR_STAR.LocalMatrix(),
          (R)1,  s21_MC_STAR.LocalMatrix() );
        // Sum the partial contributions into x21
        x21.SumScatterFrom( s21_MC_STAR );
        Scal( tauP, x21.LocalMatrix() );
        //--------------------------------------------------------------------//
        // Auxilliary alignments
        sB1_MR_STAR.FreeAlignments();
        s21_MC_STAR.FreeAlignments();
        s01_MR_STAR.FreeAlignments();
        s01_MC_STAR.FreeAlignments();
        q21_MR_STAR.FreeAlignments();
        q21_MC_STAR.FreeAlignments();
        q21_MR_MC.FreeAlignments();
        q21.FreeAlignments();
        z21_MR_STAR.FreeAlignments();
        z21_MC_STAR.FreeAlignments();
        z01_MR_STAR.FreeAlignments();
        z01_MC_STAR.FreeAlignments();
        uB1_MC_STAR.FreeAlignments();

        // Main alignments
        y10_STAR_MR.FreeAlignments();
        x10_STAR_MC.FreeAlignments();
        a12_STAR_MC.FreeAlignments();
        a10_STAR_MR.FreeAlignments();
        a01_MR_STAR.FreeAlignments();

        SlidePartitionDown
        ( eT,  e0,
               epsilon1,
         /**/ /********/
          eB,  e2 );

        SlidePartitionDown
        ( dT,  d0,
               delta1,
         /**/ /******/
          dB,  d2 );

        SlidePartitionDownDiagonal
        ( YTL, /**/ YTR,  Y00, y01,   /**/ Y02,
               /**/       y10, psi11, /**/ y12,
         /*************/ /********************/
          YBL, /**/ YBR,  Y20, y21,   /**/ Y22 );

        SlidePartitionDownDiagonal
        ( XTL, /**/ XTR,  X00, x01,   /**/ X02,
               /**/       x10, chi11, /**/ x12,
         /*************/ /********************/
          XBL, /**/ XBR,  X20, x21,   /**/ X22 );

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
    }
    PopBlocksizeStack();

    // Put back d and e
    ATL.SetDiagonal( d, 0 );
    DistMatrix<R,MC,MR> ATLExpanded(g);
    ATLExpanded.View( A, 0, 0, ATL.Height(), ATL.Width()+1 );
    ATLExpanded.SetDiagonal( e, 1 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
PanelBidiagU
( DistMatrix<Complex<R>,MC,  MR  >& A, 
  DistMatrix<Complex<R>,MD,  STAR>& tP,
  DistMatrix<Complex<R>,MD,  STAR>& tQ,
  DistMatrix<Complex<R>,MC,  MR  >& X, 
  DistMatrix<Complex<R>,MC,  MR  >& Y,
  DistMatrix<Complex<R>,MC,  STAR>& AColPan_MC_STAR,
  DistMatrix<Complex<R>,STAR,MR  >& ARowPan_STAR_MR )
{
    const int panelSize = X.Width();
#ifndef RELEASE
    PushCallStack("internal::BidiagU");
    if( A.Grid() != tP.Grid() || tP.Grid() != tQ.Grid() || 
        tQ.Grid() != X.Grid() || X.Grid() != Y.Grid() ||
        Y.Grid() != AColPan_MC_STAR.Grid() || 
        Y.Grid() != ARowPan_STAR_MR.Grid() )
        throw std::logic_error("Grids must match");
    if( tP.Height() != panelSize || tP.Width() != 1 )
        throw std::logic_error("tP was not the right size");
    if( tQ.Height() != panelSize || tQ.Width() != 1 )
        throw std::logic_error("tQ was not the right size");
    if( A.Height() < A.Width() )
        throw std::logic_error("A must be at least as tall as it is wide");
    if( A.Height() != X.Height() )
        throw std::logic_error("A and X must be the same height");
    if( A.Width() != Y.Height() )
        throw std::logic_error("Y must be the same height as A's width");
    if( X.Height() < panelSize )
        throw std::logic_error("X must be a column panel");
    if( Y.Width() != panelSize )
        throw std::logic_error("Y is the wrong width");
    if( A.ColAlignment() != X.ColAlignment() || 
        A.RowAlignment() != X.RowAlignment() )
        throw std::logic_error("A and X must be aligned");
    if( A.ColAlignment() != Y.ColAlignment() ||
        A.RowAlignment() != Y.RowAlignment() )
        throw std::logic_error("A and Y must be aligned");
#endif
    typedef Complex<R> C;

    const Grid& g = A.Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int p = g.Size();

    // Matrix views 
    DistMatrix<C,MC,MR> 
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  aB1(g), AB2(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),  alpha12L(g), a12R(g),
                         A20(g), a21(g),     A22(g),  A2L(g);
    DistMatrix<C,MC,MR>
        XTL(g), XTR(g),  X00(g), x01(g),   X02(g), 
        XBL(g), XBR(g),  x10(g), chi11(g), x12(g), 
                         X20(g), x21(g),   X22(g);
    DistMatrix<C,MC,MR>
        YTL(g), YTR(g),  Y00(g), y01(g),   Y02(g),
        YBL(g), YBR(g),  y10(g), psi11(g), y12(g),
                         Y20(g), y21(g),   Y22(g),  Y2L(g);
    DistMatrix<R,MD,STAR> d(g), dT(g), d0(g),
                                dB(g), delta1(g),
                                       d2(g);
    DistMatrix<R,MD,STAR> e(g), eT(g), e0(g),
                                eB(g), epsilon1(g),
                                       e2(g);
    DistMatrix<C,MD,STAR> tPT(g), tP0(g),
                          tPB(g), tauP1(g),
                                  tP2(g);
    DistMatrix<C,MD,STAR> tQT(g), tQ0(g),
                          tQB(g), tauQ1(g),
                                  tQ2(g);
    DistMatrix<C,MC,STAR> aB1_MC_STAR(g);
    DistMatrix<C,STAR,MR> a12_STAR_MR(g);

    // Temporary distributions
    DistMatrix<C,MR,  STAR> a01_MR_STAR(g);
    DistMatrix<C,STAR,MR  > a10_STAR_MR(g);
    DistMatrix<C,STAR,MC  > a12_STAR_MC(g);
    DistMatrix<C,STAR,MC  > x10_STAR_MC(g);
    DistMatrix<C,STAR,MR  > y10_STAR_MR(g);
    DistMatrix<C,MC,  STAR> uB1_MC_STAR(g);
    DistMatrix<C,MR,  MC  > z01_MR_MC(g);
    DistMatrix<C,MC,  STAR> z01_MC_STAR(g);
    DistMatrix<C,MR,  STAR> z01_MR_STAR(g);
    DistMatrix<C,MR,  MC  > z21_MR_MC(g);
    DistMatrix<C,MC,  STAR> z21_MC_STAR(g);
    DistMatrix<C,MR,  STAR> z21_MR_STAR(g);
    DistMatrix<C,MC,  MR  > q21(g);
    DistMatrix<C,MR,  MC  > q21_MR_MC(g);
    DistMatrix<C,MC,  STAR> q21_MC_STAR(g);
    DistMatrix<C,MR,  STAR> q21_MR_STAR(g);
    DistMatrix<C,MC,  MR  > s01(g);
    DistMatrix<C,MC,  STAR> s01_MC_STAR(g);
    DistMatrix<C,MR,  STAR> s01_MR_STAR(g);
    DistMatrix<C,MC,  STAR> s21_MC_STAR(g);
    DistMatrix<C,MR,  STAR> sB1_MR_STAR(g);

    d.AlignWithDiagonal( A, 0 );
    e.AlignWithDiagonal( A, 1 );
    d.ResizeTo( panelSize, 1 );
    e.ResizeTo( panelSize, 1 );

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDownDiagonal
    ( X, XTL, XTR,
         XBL, XBR, 0 );
    PartitionDownDiagonal
    ( Y, YTL, YTR,
         YBL, YBR, 0 );
    PartitionDown
    ( d, dT,
         dB, 0 );
    PartitionDown
    ( e, eT,
         eB, 0 );
    PartitionDown
    ( tP, tPT,
          tPB, 0 );
    PartitionDown
    ( tQ, tQT,
          tQB, 0 );
    PushBlocksizeStack( 1 );
    while( ATL.Width() < panelSize )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        RepartitionDownDiagonal
        ( XTL, /**/ XTR,  X00, /**/ x01,   X02,
         /*************/ /********************/
               /**/       x10, /**/ chi11, x12,
          XBL, /**/ XBR,  X20, /**/ x21,   X22 );
        
        RepartitionDownDiagonal
        ( YTL, /**/ YTR,  Y00, /**/ y01,   Y02,
         /*************/ /********************/
               /**/       y10, /**/ psi11, y12,
          YBL, /**/ YBR,  Y20, /**/ y21,   Y22 );

        RepartitionDown
        ( dT,  d0,
         /**/ /******/
               delta1,
          dB,  d2 );

        RepartitionDown
        ( eT,  e0,
         /**/ /********/
               epsilon1,
          eB,  e2 );

        RepartitionDown
        ( tPT,  tP0,
         /***/ /*****/
                tauP1,
          tPB,  tP2 );

        RepartitionDown
        ( tQT,  tQ0,
         /***/ /*****/
                tauQ1,
          tQB,  tQ2 );

        PartitionRight( ABR, aB1, AB2, 1 );
        PartitionRight( a12, alpha12L, a12R, 1 );

        A2L.View1x2( A20, a21 );
        Y2L.View1x2( Y20, y21 );

        a12_STAR_MR.View
        ( ARowPan_STAR_MR, ATL.Height(), ATL.Width()+1, 1, a12.Width() );
        aB1_MC_STAR.View
        ( AColPan_MC_STAR, ATL.Height(), ATL.Width(), ABR.Height(), 1 );

        // Main alignments
        a01_MR_STAR.AlignWith( ABL );
        a10_STAR_MR.AlignWith( Y20 );
        a12_STAR_MC.AlignWith( Y2L );
        x10_STAR_MC.AlignWith( A02 );
        y10_STAR_MR.AlignWith( ABL );

        // Auxilliary alignments
        uB1_MC_STAR.AlignWith( ABL );
        z01_MC_STAR.AlignWith( A02 );
        z01_MR_STAR.AlignWith( ABL );
        z21_MC_STAR.AlignWith( Y20 );
        z21_MR_STAR.AlignWith( AB2 );
        q21.AlignWith( y21 );
        q21_MR_MC.AlignWith( a12 );
        q21_MC_STAR.AlignWith( Y20 );
        q21_MR_STAR.AlignWith( A02 );
        s01_MC_STAR.AlignWith( A02 );
        s01_MR_STAR.AlignWith( X20 );
        s21_MC_STAR.AlignWith( A22 );
        sB1_MR_STAR.AlignWith( Y2L );
        
        // Auxilliary resizes
        uB1_MC_STAR.ResizeTo( ABL.Height(), 1 );
        z01_MR_STAR.ResizeTo( A00.Width(), 1 );
        z21_MC_STAR.ResizeTo( A22.Width(), 1 );
        z21_MR_STAR.ResizeTo( A22.Width(), 1 );
        q21_MC_STAR.ResizeTo( A22.Width(), 1 );
        q21_MR_STAR.ResizeTo( A22.Width(), 1 );
        s01_MC_STAR.ResizeTo( A00.Height(), 1 );
        s21_MC_STAR.ResizeTo( A22.Height(), 1 );
        sB1_MR_STAR.ResizeTo( Y2L.Width(), 1 );

        const bool thisIsMyRow = ( g.Row() == alpha11.ColAlignment() );
        const bool thisIsMyCol = ( g.Col() == alpha11.RowAlignment() );
        const bool nextIsMyCol = ( g.Col() == a12.RowAlignment() ) ;
        const bool firstIteration = ( ATL.Height() == 0 );
        //--------------------------------------------------------------------//

        // Update the current column of A:
        //   aB1 := aB1 - ABL y10^H - XBL a01
        if( !firstIteration )
        {
            Conjugate( y10 );
            y10_STAR_MR = y10;
            // uB1[MC,* ] := ABL[MC,MR] y10^H[MR,* ]
            a01_MR_STAR = a01;
            Gemv
            ( NORMAL, 
              (C)1, ABL.LocalMatrix(), y10_STAR_MR.LocalMatrix(), 
              (C)0, uB1_MC_STAR.LocalMatrix() );
            // uB1[MC,* ] := uB1[MC,* ] + XBL[MC,MR] a01[MR,* ]
            //             = ABL[MC,MR] y10^H[MR,* ] + XBL[MC,MR] a01[MR,* ]
            Gemv
            ( NORMAL,
              (C)1, XBL.LocalMatrix(), a01_MR_STAR.LocalMatrix(),
              (C)1, uB1_MC_STAR.LocalMatrix() );
            // Sum the partial contributions and subtract from aB1
            aB1.SumScatterUpdate( (C)-1, uB1_MC_STAR );
        }

        // Find tauQ, u, and delta such that
        //     I - conj(tauQ) | 1 | | 1, u^H | | alpha11 | = | delta |
        //                    | u |            |   a21   | = |   0   |
        C tauQ = 0;
        if( thisIsMyCol )
        {
            tauQ = ColReflector( alpha11, a21 );
            if( thisIsMyRow )
            {
                tauQ1.SetLocalEntry(0,0,tauQ);
                // Store delta and force | alpha11 | = | 1 |
                //                       |   a21   |   | u |
                delta1.SetLocalEntry(0,0,alpha11.GetRealLocalEntry(0,0));
                alpha11.SetLocalEntry(0,0,(C)1);
            }
        }

        //
        // y21 := tauQ ( AB2^H aB1 - A02^H XBL^H aB1 - Y20 ABL^H aB1 )
        //
        aB1_MC_STAR = aB1;
        // z01[MR,* ] := ABL^H[MR,MC] aB1[MC,* ]
        Gemv
        ( ADJOINT,
          (C)1, ABL.LocalMatrix(), aB1_MC_STAR.LocalMatrix(),
          (C)0, z01_MR_STAR.LocalMatrix() );
        // z21[MR,* ] := AB2^H[MR,MC] aB1[MC,* ]
        Gemv
        ( ADJOINT,
          (C)1, AB2.LocalMatrix(), aB1_MC_STAR.LocalMatrix(),
          (C)0, z21_MR_STAR.LocalMatrix() );
        // Sum the partial contributions
        z01_MR_STAR.SumOverCol();
        // z21[MC,* ] := Y20[MC,MR] z01[MR,* ] = Y20[MC,MR] (ABL^H aB1)[MR,* ]
        Gemv
        ( NORMAL,
          (C)1, Y20.LocalMatrix(), z01_MR_STAR.LocalMatrix(),
          (C)0, z21_MC_STAR.LocalMatrix() );
        // z01[MR,* ] := XBL^H[MR,MC] aB1[MC,* ]
        Gemv
        ( ADJOINT,
          (C)1, XBL.LocalMatrix(), aB1_MC_STAR.LocalMatrix(),
          (C)0, z01_MR_STAR.LocalMatrix() );
        // Sum the partial contributions to z01[MR,* ] and scatter the result
        z01_MR_MC.SumScatterFrom( z01_MR_STAR );
        // Redistribute the scattered summation
        z01_MC_STAR = z01_MR_MC;
        // z21[MR,* ] := z21[MR,* ] - A02^H[MR,MC] z01[MC,* ]
        //             = AB2^H[MR,MC] aB1[MC,* ] - 
        //               A02^H[MR,MC] (XBL^H aB1)[MC,* ]
        Gemv
        ( ADJOINT,
          (C)-1, A02.LocalMatrix(), z01_MC_STAR.LocalMatrix(),
          (C)1, z21_MR_STAR.LocalMatrix() );
        // Sum the partial contributions to z21[MR,* ] and scatter the result
        z21_MR_MC.SumScatterFrom( z21_MR_STAR );
        // Redistribute (and rename) the scattered summation
        y21 = z21_MR_MC;
        // Subtract z21 = Y20 ABL^H aB1 from y21
        y21.SumScatterUpdate( (C)-1, z21_MC_STAR );
        if( thisIsMyCol )
            Scal( tauQ, y21 );

        //
        // y21 := y21 + Y20 a10^H
        //
        Conjugate( a10 );
        a10_STAR_MR = a10;
        Conjugate( a10 );
        // q21[MC,* ] := Y20[MC,MR] a10^H[MR,* ]
        Gemv
        ( NORMAL, 
          (C)1, Y20.LocalMatrix(), a10_STAR_MR.LocalMatrix(),
          (C)0, q21_MC_STAR.LocalMatrix() );
        // Sum the partial contributions
        q21.SumScatterFrom( q21_MC_STAR );
        if( thisIsMyCol )
            Axpy( (C)1, y21, q21 );

        //
        // a12 := conj(a12 - a10 Y20^H - x10 A02)
        //
        Conjugate( x10 );
        x10_STAR_MC = x10;
        Conjugate( x10 );
        q21_MR_MC = q21;
        // q21[MR,* ] := A02^H[MR,MC] x10^H[MC,* ]
        Gemv
        ( ADJOINT,
          (C)1, A02.LocalMatrix(), x10_STAR_MC.LocalMatrix(),
          (C)0, q21_MR_STAR.LocalMatrix() );
        // Sum the partial contributions onto q21[MR,MC] = (Y20 a10^H)[MR,MC]
        q21_MR_MC.SumScatterUpdate( (C)1, q21_MR_STAR );
        // a12 := conj(a12) - q21^T = conj(a12 - a10 Y20^H - x10 A02)
        Conjugate( a12 );
        if( thisIsMyRow )
        {
            const int localWidth = a12.LocalWidth();
            C* a12Buffer = a12.LocalBuffer();
            const C* q21Buffer = q21_MR_MC.LockedLocalBuffer();
            const int a12LDim = a12.LocalLDim();
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                a12Buffer[jLocal*a12LDim] -= q21Buffer[jLocal];
        }

        // Find tauP, v, and epsilon such that
        //     I - conj(tauP) | 1 | | 1, v^H | | alpha12L | = | epsilon |
        //                    | v |            |  a12R^T  |   |    0    |
        C tauP = 0;
        if( thisIsMyRow )
        {
            tauP = RowReflector( alpha12L, a12R );
            if( nextIsMyCol )
            {
                tauP1.SetLocalEntry(0,0,tauP);
                // Store epsilon and force | alpha12L | = | 1 |
                //                         |  a12R^T  |   | v |
                epsilon1.SetLocalEntry(0,0,alpha12L.GetRealLocalEntry(0,0));
                alpha12L.SetLocalEntry(0,0,(C)1);
            }
        }
        mpi::Broadcast( &tauP, 1, alpha11.ColAlignment(), g.ColComm() );

        //
        // (Keep in mind that a12 is currently overwritten with its conjugate.
        //  We will use the 'true' value in the following comments.)
        //
        // x21 := conj(tauP) ( A22 a12^H - A2L Y2L^H a12^H - X20 A02 a12^H )
        //
        a12_STAR_MR = a12;
        a12_STAR_MC = a12;
        // s21[MC,* ] := A22[MC,MR] a12^H[MR,* ]
        Gemv
        ( NORMAL,
          (C)1, A22.LocalMatrix(), a12_STAR_MR.LocalMatrix(),
          (C)0, s21_MC_STAR.LocalMatrix() );
        // sB1[MR,* ] := Y2L^H[MR,MC] a12^H[MC,* ]
        Gemv
        ( ADJOINT,
          (C)1, Y2L.LocalMatrix(), a12_STAR_MC.LocalMatrix(),
          (C)0, sB1_MR_STAR.LocalMatrix() );
        // Sum the partial contributions
        sB1_MR_STAR.SumOverCol(); 
        // s21[MC,* ] := s21[MC,* ] - A2L[MC,MR] sB1[MR,* ]
        //             = A22[MC,MR] a12^H[MR,* ] - A2L[MC,MR] sB1[MR,* ]
        // (still needs to be summed within each process row)
        Gemv
        ( NORMAL, 
          (C)-1, A2L.LocalMatrix(), sB1_MR_STAR.LocalMatrix(),
          (C)1,  s21_MC_STAR.LocalMatrix() );
        // s01[MC,* ] := A02[MC,MR] a12^H[MR,* ]
        Gemv
        ( NORMAL,
          (C)1, A02.LocalMatrix(), a12_STAR_MR.LocalMatrix(),
          (C)0, s01_MC_STAR.LocalMatrix() );
        // Sum the partial contributions and then redistribute
        s01.SumScatterFrom( s01_MC_STAR ); // TODO: SumScatter to [VC,* ]?
        s01_MR_STAR = s01;
        // s21[MC,* ] := s21[MC,* ] - X20[MC,MR] s01[MR,* ]
        //             = A22[MC,MR] a12^H[MR,* ] - A2L[MC,MR] sB1[MR,* ]
        //                                       - X20[MC,MR] s01[MR,* ]
        Gemv
        ( NORMAL,
          (C)-1, X20.LocalMatrix(), s01_MR_STAR.LocalMatrix(),
          (C)1,  s21_MC_STAR.LocalMatrix() );
        // Sum the partial contributions into x21
        x21.SumScatterFrom( s21_MC_STAR );
        Scal( tauP, x21.LocalMatrix() );

        // Undo the in-place conjugation of a12
        Conjugate( a12 );
        Conjugate( a12_STAR_MR );
        //--------------------------------------------------------------------//
        // Auxilliary alignments
        sB1_MR_STAR.FreeAlignments();
        s21_MC_STAR.FreeAlignments();
        s01_MR_STAR.FreeAlignments();
        s01_MC_STAR.FreeAlignments();
        q21_MR_STAR.FreeAlignments();
        q21_MC_STAR.FreeAlignments();
        q21_MR_MC.FreeAlignments();
        q21.FreeAlignments();
        z21_MR_STAR.FreeAlignments();
        z21_MC_STAR.FreeAlignments();
        z01_MR_STAR.FreeAlignments();
        z01_MC_STAR.FreeAlignments();
        uB1_MC_STAR.FreeAlignments();

        // Main alignments
        y10_STAR_MR.FreeAlignments();
        x10_STAR_MC.FreeAlignments();
        a12_STAR_MC.FreeAlignments();
        a10_STAR_MR.FreeAlignments();
        a01_MR_STAR.FreeAlignments();

        SlidePartitionDown
        ( tQT,  tQ0,
                tauQ1,
         /***/ /*****/
          tQB,  tQ2 );

        SlidePartitionDown
        ( tPT,  tP0,
                tauP1,
         /***/ /*****/
          tPB,  tP2 );

        SlidePartitionDown
        ( eT,  e0,
               epsilon1,
         /**/ /********/
          eB,  e2 );

        SlidePartitionDown
        ( dT,  d0,
               delta1,
         /**/ /******/
          dB,  d2 );

        SlidePartitionDownDiagonal
        ( YTL, /**/ YTR,  Y00, y01,   /**/ Y02,
               /**/       y10, psi11, /**/ y12,
         /*************/ /********************/
          YBL, /**/ YBR,  Y20, y21,   /**/ Y22 );

        SlidePartitionDownDiagonal
        ( XTL, /**/ XTR,  X00, x01,   /**/ X02,
               /**/       x10, chi11, /**/ x12,
         /*************/ /********************/
          XBL, /**/ XBR,  X20, x21,   /**/ X22 );

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
    }
    PopBlocksizeStack();

    // Put back d and e
    ATL.SetRealDiagonal( d, 0 );
    DistMatrix<Complex<R>,MC,MR> ATLExpanded(g);
    ATLExpanded.View( A, 0, 0, ATL.Height(), ATL.Width()+1 );
    ATLExpanded.SetRealDiagonal( e, 1 );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
