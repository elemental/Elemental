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

template<typename R> // representation of a real number 
inline void
elemental::advanced::internal::HermitianTridiagLSquare
( DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::HermitianTridiagLSquare");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Grid().Height() != A.Grid().Width() )
        throw std::logic_error("The process grid must be square");
#endif
    const Grid& g = A.Grid();

    if( g.InGrid() )
    {
        // Matrix views 
        DistMatrix<R,MC,MR> 
            ATL(g), ATR(g),  A00(g), A01(g), A02(g), 
            ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                             A20(g), A21(g), A22(g);

        // Temporary distributions
        DistMatrix<R,STAR,STAR> A11_STAR_STAR(g);
        DistMatrix<R,MC,  STAR> APan_MC_STAR(g);
        DistMatrix<R,MR,  STAR> APan_MR_STAR(g);
        DistMatrix<R,MC,  STAR> A11_MC_STAR(g);
        DistMatrix<R,MR,  STAR> A11_MR_STAR(g);
        DistMatrix<R,MC,  STAR> A21_MC_STAR(g);
        DistMatrix<R,MR,  STAR> A21_MR_STAR(g);
        DistMatrix<R,MC,  MR  > WPan(g);
        DistMatrix<R,MC,  STAR> WPan_MC_STAR(g);
        DistMatrix<R,MR,  STAR> WPan_MR_STAR(g);
        DistMatrix<R,MC,  STAR> W11_MC_STAR(g);
        DistMatrix<R,MR,  STAR> W11_MR_STAR(g);
        DistMatrix<R,MC,  STAR> W21_MC_STAR(g);
        DistMatrix<R,MR,  STAR> W21_MR_STAR(g);

        PartitionDownDiagonal
        ( A, ATL, ATR,
             ABL, ABR, 0 );
        while( ATL.Height() < A.Height() )
        {
            RepartitionDownDiagonal
            ( ATL, /**/ ATR,  A00, /**/ A01, A02,
             /*************/ /******************/
                   /**/       A10, /**/ A11, A12,
              ABL, /**/ ABR,  A20, /**/ A21, A22 );

            if( A22.Height() > 0 )
            {
                APan_MC_STAR.AlignWith( A11 );
                APan_MR_STAR.AlignWith( A11 );
                APan_MC_STAR.ResizeTo( ABR.Height(), A11.Width() );
                APan_MR_STAR.ResizeTo( ABR.Height(), A11.Width() );
                WPan.AlignWith( A11 );
                WPan_MC_STAR.AlignWith( A11 );
                WPan_MR_STAR.AlignWith( A11 );
                WPan.ResizeTo( ABR.Height(), A11.Width() );
                WPan_MC_STAR.ResizeTo( ABR.Height(), A11.Width() );
                WPan_MR_STAR.ResizeTo( ABR.Height(), A11.Width() );
                PartitionDown
                ( APan_MC_STAR, A11_MC_STAR,
                                A21_MC_STAR, A11.Height() );
                PartitionDown
                ( APan_MR_STAR, A11_MR_STAR,
                                A21_MR_STAR, A11.Height() );
                PartitionDown
                ( WPan_MC_STAR, W11_MC_STAR,
                                W21_MC_STAR, A11.Height() );
                PartitionDown
                ( WPan_MR_STAR, W11_MR_STAR,
                                W21_MR_STAR, A11.Height() );
                //------------------------------------------------------------//
                // Accumulate the Householder vectors into A21 and form W21 
                // such that subtracting (A21 W21' + W21 A21') is equal to 
                // successively applying the similarity transformations 
                // (I-tau h h')A22(I-tau h h') for each (tau,h).
                //
                // APan[MC,* ], APan[MR,* ], WPan[MC,* ], and WPan[MR,* ] are 
                // formed during the panel factorization.
                advanced::internal::HermitianPanelTridiagLSquare
                ( ABR, WPan, 
                  APan_MC_STAR, APan_MR_STAR, WPan_MC_STAR, WPan_MR_STAR );
                basic::internal::LocalTrr2k
                ( LOWER, TRANSPOSE, TRANSPOSE,
                  (R)-1, A21_MC_STAR, W21_MR_STAR,
                         W21_MC_STAR, A21_MR_STAR,
                  (R)1, A22 );
                //------------------------------------------------------------//
                APan_MC_STAR.FreeAlignments();
                APan_MR_STAR.FreeAlignments();
                WPan.FreeAlignments();
                WPan_MC_STAR.FreeAlignments();
                WPan_MR_STAR.FreeAlignments();
            }
            else
            {
                A11_STAR_STAR = A11;
                advanced::HermitianTridiag
                ( LOWER, A11_STAR_STAR.LocalMatrix() );
                A11 = A11_STAR_STAR;
            }

            SlidePartitionDownDiagonal
            ( ATL, /**/ ATR,  A00, A01, /**/ A02,
                   /**/       A10, A11, /**/ A12,
             /*************/ /******************/
              ABL, /**/ ABR,  A20, A21, /**/ A22 );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
inline void
elemental::advanced::internal::HermitianTridiagLSquare
( DistMatrix<std::complex<R>,MC,  MR  >& A,
  DistMatrix<std::complex<R>,STAR,STAR>& t )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::HermitianTridiagLSquare");
    if( A.Grid() != t.Grid() )
        throw std::logic_error("{A,t} must be distributed over the same grid");
#endif
    const Grid& g = A.Grid();
#ifndef RELEASE
    if( g.Height() != g.Width() )
        throw std::logic_error("The process grid must be square");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( t.Viewing() )
        throw std::logic_error("t must not be a view");
#endif
    typedef std::complex<R> C;

    DistMatrix<C,MD,STAR> tDiag(g);
    tDiag.AlignWithDiagonal( A, -1 );
    tDiag.ResizeTo( A.Height()-1, 1 );

    if( g.InGrid() )
    {
        // Matrix views 
        DistMatrix<C,MC,MR> 
            ATL(g), ATR(g),  A00(g), A01(g), A02(g), 
            ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                             A20(g), A21(g), A22(g);
        DistMatrix<C,MD,STAR> tT(g),  t0(g), 
                              tB(g),  t1(g),
                                      t2(g);

        // Temporary distributions
        DistMatrix<C,STAR,STAR> A11_STAR_STAR(g);
        DistMatrix<C,MC,  STAR> APan_MC_STAR(g);
        DistMatrix<C,MR,  STAR> APan_MR_STAR(g);
        DistMatrix<C,MC,  STAR> A11_MC_STAR(g);
        DistMatrix<C,MR,  STAR> A11_MR_STAR(g);
        DistMatrix<C,MC,  STAR> A21_MC_STAR(g);
        DistMatrix<C,MR,  STAR> A21_MR_STAR(g);
        DistMatrix<C,MC,  MR  > WPan(g);
        DistMatrix<C,MC,  STAR> WPan_MC_STAR(g);
        DistMatrix<C,MR,  STAR> WPan_MR_STAR(g);
        DistMatrix<C,MC,  STAR> W11_MC_STAR(g);
        DistMatrix<C,MR,  STAR> W11_MR_STAR(g);
        DistMatrix<C,MC,  STAR> W21_MC_STAR(g);
        DistMatrix<C,MR,  STAR> W21_MR_STAR(g);
        DistMatrix<C,STAR,STAR> t1_STAR_STAR(g);

        PartitionDownDiagonal
        ( A, ATL, ATR,
             ABL, ABR, 0 );
        PartitionDown
        ( tDiag, tT,
                 tB, 0 );
        while( ATL.Height() < A.Height() )
        {
            RepartitionDownDiagonal
            ( ATL, /**/ ATR,  A00, /**/ A01, A02,
             /*************/ /******************/
                   /**/       A10, /**/ A11, A12,
              ABL, /**/ ABR,  A20, /**/ A21, A22 );

            RepartitionDown
            ( tT,  t0,
             /**/ /**/
                   t1,
              tB,  t2 );
            
            if( A22.Height() > 0 )
            {
                APan_MC_STAR.AlignWith( A11 );
                APan_MR_STAR.AlignWith( A11 );
                APan_MC_STAR.ResizeTo( ABR.Height(), A11.Width() );
                APan_MR_STAR.ResizeTo( ABR.Height(), A11.Width() );
                WPan.AlignWith( A11 );
                WPan_MC_STAR.AlignWith( A11 );
                WPan_MR_STAR.AlignWith( A11 );
                WPan.ResizeTo( ABR.Height(), A11.Width() );
                WPan_MC_STAR.ResizeTo( ABR.Height(), A11.Width() );
                WPan_MR_STAR.ResizeTo( ABR.Height(), A11.Width() );
                PartitionDown
                ( APan_MC_STAR, A11_MC_STAR,
                                A21_MC_STAR, A11.Height() );
                PartitionDown
                ( APan_MR_STAR, A11_MR_STAR,
                                A21_MR_STAR, A11.Height() );
                PartitionDown
                ( WPan_MC_STAR, W11_MC_STAR,
                                W21_MC_STAR, A11.Height() );
                PartitionDown
                ( WPan_MR_STAR, W11_MR_STAR,
                                W21_MR_STAR, A11.Height() );
                //------------------------------------------------------------//
                // Accumulate the Householder vectors into A21 and form W21 
                // such that subtracting (A21 W21' + W21 A21') is equal to 
                // successively applying the similarity transformations 
                // (I-conj(tau) h h')A22(I-tau h h') for each (tau,h).
                //
                // APan[MC,* ], APan[MR,* ], WPan[MC,* ], and WPan[MR,* ] are 
                // formed during the panel factorization.
                advanced::internal::HermitianPanelTridiagLSquare
                ( ABR, WPan, t1,
                  APan_MC_STAR, APan_MR_STAR, WPan_MC_STAR, WPan_MR_STAR );
                basic::internal::LocalTrr2k
                ( LOWER, ADJOINT, ADJOINT,
                  (C)-1, A21_MC_STAR, W21_MR_STAR,
                         W21_MC_STAR, A21_MR_STAR,
                  (C)1, A22 );
                //------------------------------------------------------------//
                APan_MC_STAR.FreeAlignments();
                APan_MR_STAR.FreeAlignments();
                WPan.FreeAlignments();
                WPan_MC_STAR.FreeAlignments();
                WPan_MR_STAR.FreeAlignments();
            }
            else
            {
                A11_STAR_STAR = A11;
                t1_STAR_STAR.ResizeTo( t1.Height(), 1 );

                advanced::HermitianTridiag
                ( LOWER, A11_STAR_STAR.LocalMatrix(), 
                  t1_STAR_STAR.LocalMatrix() );

                A11 = A11_STAR_STAR;
                t1 = t1_STAR_STAR;
            }

            SlidePartitionDownDiagonal
            ( ATL, /**/ ATR,  A00, A01, /**/ A02,
                   /**/       A10, A11, /**/ A12,
             /*************/ /******************/
              ABL, /**/ ABR,  A20, A21, /**/ A22 );

            SlidePartitionDown
            ( tT,  t0,
                   t1,
             /**/ /**/
              tB,  t2 );
        }
    }
    // Redistribute from matrix-diagonal form to fully replicated
    t = tDiag;
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX
