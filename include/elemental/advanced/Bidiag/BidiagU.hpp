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

template<typename R>
inline void 
elemental::advanced::internal::BidiagU( DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::BidiagU");
    if( A.Height() < A.Width() )
        throw std::logic_error("A must be at least as tall as it is wide");
#endif
    const Grid& g = A.Grid();

    // Matrix views 
    DistMatrix<R,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g), 
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary distributions
    DistMatrix<R,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<R,MC,  MR  > XColPan(g);
    DistMatrix<R,MC,  MR  > YRowPan(g);
    DistMatrix<R,MC,  STAR> AColPan_MC_STAR(g), A11_MC_STAR(g),
                                                A21_MC_STAR(g);
    DistMatrix<R,MC,  STAR> XColPan_MC_STAR(g), X11_MC_STAR(g),
                                                X21_MC_STAR(g);
    DistMatrix<R,STAR,MR  > ARowPan_STAR_MR(g), A11_STAR_MR(g), A12_STAR_MR(g);
    DistMatrix<R,STAR,MR  > YRowPan_STAR_MR(g), Y11_STAR_MR(g), Y12_STAR_MR(g);

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        if( A22.Width() > 0 )
        {
            XColPan.AlignWith( A11 );
            YRowPan.AlignWith( A11 );
            AColPan_MC_STAR.AlignWith( A11 );
            XColPan_MC_STAR.AlignWith( A11 );
            ARowPan_STAR_MR.AlignWith( A11 );
            YRowPan_STAR_MR.AlignWith( A11 );
            //----------------------------------------------------------------//
            XColPan.ResizeTo( ABR.Height(), A11.Width() );
            YRowPan.ResizeTo( A11.Height(), ABR.Width() );
            AColPan_MC_STAR.ResizeTo( ABR.Height(), A11.Width() );
            XColPan_MC_STAR.ResizeTo( ABR.Height(), A11.Width() );
            ARowPan_STAR_MR.ResizeTo( A11.Height(), ABR.Width() );
            YRowPan_STAR_MR.ResizeTo( A11.Height(), ABR.Width() );

            advanced::internal::PanelBidiagU
            ( ABR, XColPan, YRowPan, 
              AColPan_MC_STAR, XColPan_MC_STAR,
              ARowPan_STAR_MR, YRowPan_STAR_MR );

            PartitionDown
            ( AColPan_MC_STAR, A11_MC_STAR,
                               A21_MC_STAR, A11.Height() );
            PartitionDown
            ( XColPan_MC_STAR, X11_MC_STAR,
                               X21_MC_STAR, A11.Height() );
            PartitionRight
            ( ARowPan_STAR_MR, A11_STAR_MR, A12_STAR_MR, A11.Width() );
            PartitionRight
            ( YRowPan_STAR_MR, Y11_STAR_MR, Y12_STAR_MR, A11.Width() );

            basic::internal::LocalGemm
            ( NORMAL, NORMAL, (R)-1, A21_MC_STAR, YRowPan_STAR_MR, (R)1, A22 );
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, (R)-1, XColPan_MC_STAR, A12_STAR_MR, (R)1, A22 );
            //----------------------------------------------------------------//
            YRowPan_STAR_MR.FreeAlignments();
            ARowPan_STAR_MR.FreeAlignments();
            XColPan_MC_STAR.FreeAlignments();
            AColPan_MC_STAR.FreeAlignments();
            YRowPan.FreeAlignments();
            XColPan.FreeAlignments();
        }
        else
        {
            A11_STAR_STAR = A11;
            advanced::Bidiag( A11_STAR_STAR.LocalMatrix() );
            A11 = A11_STAR_STAR;
        }

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
elemental::advanced::internal::BidiagU
( DistMatrix<std::complex<R>,MC,  MR  >& A,
  DistMatrix<std::complex<R>,STAR,STAR>& tP,
  DistMatrix<std::complex<R>,STAR,STAR>& tQ )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::BidiagU");
    if( A.Grid() != tP.Grid() || tP.Grid() != tQ.Grid() )
        throw std::logic_error
        ("{A,tP,tQ} must be distributed over the same grid");
    if( A.Height() < A.Width() )
        throw std::logic_error("A must be at least as tall as it is wide");
    if( tP.Viewing() || tQ.Viewing() )
        throw std::logic_error("tP and tQ must not be views");
#endif
    typedef std::complex<R> C;

    const Grid& g = A.Grid();
    const int tPHeight = std::max(A.Width()-1,0);
    const int tQHeight = A.Width();
    DistMatrix<C,MD,STAR> tPDiag(g), tQDiag(g);
    tPDiag.AlignWithDiagonal( A, 1 );
    tQDiag.AlignWithDiagonal( A, 0 );
    tPDiag.ResizeTo( tPHeight, 1 );
    tQDiag.ResizeTo( tQHeight, 1 );

    // Matrix views 
    DistMatrix<C,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g), 
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<C,MD,STAR> tPT(g),  tP0(g), 
                          tPB(g),  tP1(g),
                                   tP2(g);
    DistMatrix<C,MD,STAR> tQT(g),  tQ0(g), 
                          tQB(g),  tQ1(g),
                                   tQ2(g);

    // Temporary distributions
    DistMatrix<C,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<C,STAR,STAR> tP1_STAR_STAR(g);
    DistMatrix<C,STAR,STAR> tQ1_STAR_STAR(g);
    DistMatrix<C,MC,  MR  > XColPan(g);
    DistMatrix<C,MC,  MR  > YRowPan(g);
    DistMatrix<C,MC,  STAR> AColPan_MC_STAR(g), A11_MC_STAR(g),
                                                A21_MC_STAR(g);
    DistMatrix<C,MC,  STAR> XColPan_MC_STAR(g), X11_MC_STAR(g),
                                                X21_MC_STAR(g);
    DistMatrix<C,STAR,MR  > ARowPan_STAR_MR(g), A11_STAR_MR(g), A12_STAR_MR(g);
    DistMatrix<C,STAR,MR  > YRowPan_STAR_MR(g), Y11_STAR_MR(g), Y12_STAR_MR(g);

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( tPDiag, tPT,
              tPB, 0 );
    PartitionDown
    ( tQDiag, tQT,
              tQB, 0 );
    while( ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        RepartitionDown
        ( tPT,  tP0,
         /***/ /***/
                tP1,
          tPB,  tP2 );
        
        RepartitionDown
        ( tQT,  tQ0,
         /***/ /***/
                tQ1,
          tQB,  tQ2 );
        
        if( A22.Width() > 0 )
        {
            XColPan.AlignWith( A11 );
            YRowPan.AlignWith( A11 );
            AColPan_MC_STAR.AlignWith( A11 );
            XColPan_MC_STAR.AlignWith( A11 );
            ARowPan_STAR_MR.AlignWith( A11 );
            YRowPan_STAR_MR.AlignWith( A11 );
            //----------------------------------------------------------------//
            XColPan.ResizeTo( ABR.Height(), A11.Width() );
            YRowPan.ResizeTo( A11.Height(), ABR.Width() );
            AColPan_MC_STAR.ResizeTo( ABR.Height(), A11.Width() );
            XColPan_MC_STAR.ResizeTo( ABR.Height(), A11.Width() );
            ARowPan_STAR_MR.ResizeTo( A11.Height(), ABR.Width() );
            YRowPan_STAR_MR.ResizeTo( A11.Height(), ABR.Width() );

            advanced::internal::PanelBidiagU
            ( ABR, tP1, tQ1, XColPan, YRowPan, 
              AColPan_MC_STAR, XColPan_MC_STAR,
              ARowPan_STAR_MR, YRowPan_STAR_MR );

            PartitionDown
            ( AColPan_MC_STAR, A11_MC_STAR,
                               A21_MC_STAR, A11.Height() );
            PartitionDown
            ( XColPan_MC_STAR, X11_MC_STAR,
                               X21_MC_STAR, A11.Height() );
            PartitionRight
            ( ARowPan_STAR_MR, A11_STAR_MR, A12_STAR_MR, A11.Width() );
            PartitionRight
            ( YRowPan_STAR_MR, Y11_STAR_MR, Y12_STAR_MR, A11.Width() );

            basic::internal::LocalGemm
            ( NORMAL, NORMAL, (C)-1, A21_MC_STAR, YRowPan_STAR_MR, (C)1, A22 );
            basic::internal::LocalGemm
            ( NORMAL, NORMAL, (C)-1, XColPan_MC_STAR, A12_STAR_MR, (C)1, A22 );
            //----------------------------------------------------------------//
            YRowPan_STAR_MR.FreeAlignments();
            ARowPan_STAR_MR.FreeAlignments();
            XColPan_MC_STAR.FreeAlignments();
            AColPan_MC_STAR.FreeAlignments();
            YRowPan.FreeAlignments();
            XColPan.FreeAlignments();
        }
        else
        {
            A11_STAR_STAR = A11;
            tP1_STAR_STAR.ResizeTo( tP1.Height(), 1 );
            tQ1_STAR_STAR.ResizeTo( tQ1.Height(), 1 );

            advanced::Bidiag
            ( A11_STAR_STAR.LocalMatrix(), 
              tP1_STAR_STAR.LocalMatrix(), 
              tQ1_STAR_STAR.LocalMatrix() );

            A11 = A11_STAR_STAR;
            tP1 = tP1_STAR_STAR;
            tQ1 = tQ1_STAR_STAR;
        }

        SlidePartitionDown
        ( tQT,  tQ0,
                tQ1,
         /***/ /***/
          tQB,  tQ2 );

        SlidePartitionDown
        ( tPT,  tP0,
                tP1,
         /***/ /***/
          tPB,  tP2 );
        
        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }

    // Redistribute from matrix-diagonal form to fully replicated
    tP = tPDiag;
    tQ = tQDiag;
#ifndef RELEASE
    PopCallStack();
#endif
}
