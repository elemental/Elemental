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
elemental::advanced::internal::PanelQR
( DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::PanelQR");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<R,MC,MR>
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  aLeftCol(g), ARightPan(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),
                         A20(g), a21(g),     A22(g);

    // Temporary distributions
    DistMatrix<R,MC,STAR> aLeftCol_MC_STAR(g);
    DistMatrix<R,MR,STAR> Z_MR_STAR(g);

    PushBlocksizeStack( 1 );
    PartitionDownLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        aLeftCol.View2x1( alpha11,
                          a21 );

        ARightPan.View2x1( a12,
                           A22 );

        aLeftCol_MC_STAR.AlignWith( ARightPan );
        Z_MR_STAR.AlignWith( ARightPan );
        Z_MR_STAR.ResizeTo( ARightPan.Width(), 1 );
        //--------------------------------------------------------------------//
        R tau = advanced::Reflector( alpha11, a21 );

        bool myDiagonalEntry = ( g.MCRank() == alpha11.ColAlignment() && 
                                 g.MRRank() == alpha11.RowAlignment() );
        R alpha = (R)0;
        if( myDiagonalEntry )
        {
            alpha = alpha11.GetLocalEntry(0,0);
            alpha11.SetLocalEntry(0,0,1);
        }

        aLeftCol_MC_STAR = aLeftCol;

        basic::Gemv
        ( TRANSPOSE, 
          (R)1, ARightPan.LockedLocalMatrix(), 
                aLeftCol_MC_STAR.LockedLocalMatrix(),
          (R)0, Z_MR_STAR.LocalMatrix() );
        Z_MR_STAR.SumOverCol(); 

        basic::Ger
        ( -tau, 
          aLeftCol_MC_STAR.LockedLocalMatrix(), 
          Z_MR_STAR.LockedLocalMatrix(),
          ARightPan.LocalMatrix() );

        if( myDiagonalEntry )
            alpha11.SetLocalEntry(0,0,alpha);
        //--------------------------------------------------------------------//
        aLeftCol_MC_STAR.FreeAlignments();
        Z_MR_STAR.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
inline void
elemental::advanced::internal::PanelQR
( DistMatrix<std::complex<R>,MC,MR  >& A,
  DistMatrix<std::complex<R>,MD,STAR>& t )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::PanelQR");
    if( A.Grid() != t.Grid() )
        throw std::logic_error("{A,t} must be distributed over the same grid");
    if( t.Height() != min(A.Height(),A.Width()) || t.Width() != 1 )
        throw std::logic_error
        ("t must be a vector of height equal to the minimum dimension of A");
    if( !t.AlignedWithDiagonal( A, 0 ) )
        throw std::logic_error("t must be aligned with A's main diagonal");
#endif
    typedef std::complex<R> C;
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<C,MC,MR>
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  aLeftCol(g), ARightPan(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),
                         A20(g), a21(g),     A22(g);
    DistMatrix<C,MD,STAR>
        tT(g),  t0(g),
        tB(g),  tau1(g),
                t2(g);

    // Temporary distributions
    DistMatrix<C,MC,STAR> aLeftCol_MC_STAR(g);
    DistMatrix<C,MR,STAR> Z_MR_STAR(g);

    PushBlocksizeStack( 1 );
    PartitionDownLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( t, tT,
         tB, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        RepartitionDown
        ( tT,  t0,
         /**/ /****/
               tau1, 
          tB,  t2 );

        aLeftCol.View2x1( alpha11,
                          a21 );

        ARightPan.View2x1( a12,
                           A22 );

        aLeftCol_MC_STAR.AlignWith( ARightPan );
        Z_MR_STAR.AlignWith( ARightPan );
        Z_MR_STAR.ResizeTo( ARightPan.Width(), 1 );
        //--------------------------------------------------------------------//
        C tau = advanced::Reflector( alpha11, a21 );
        tau1.Set( 0, 0, tau );

        bool myDiagonalEntry = ( g.MCRank() == alpha11.ColAlignment() && 
                                 g.MRRank() == alpha11.RowAlignment() );
        C alpha = (C)0;
        if( myDiagonalEntry )
        {
            alpha = alpha11.GetLocalEntry(0,0);
            alpha11.SetLocalEntry(0,0,1);
        }

        aLeftCol_MC_STAR = aLeftCol;

        basic::Gemv
        ( ADJOINT, 
          (C)1, ARightPan.LockedLocalMatrix(), 
                aLeftCol_MC_STAR.LockedLocalMatrix(),
          (C)0, Z_MR_STAR.LocalMatrix() );
        Z_MR_STAR.SumOverCol(); 

        basic::Ger
        ( -conj(tau), 
          aLeftCol_MC_STAR.LockedLocalMatrix(), 
          Z_MR_STAR.LockedLocalMatrix(),
          ARightPan.LocalMatrix() );

        if( myDiagonalEntry )
            alpha11.SetLocalEntry(0,0,alpha);
        //--------------------------------------------------------------------//
        aLeftCol_MC_STAR.FreeAlignments();
        Z_MR_STAR.FreeAlignments();

        SlidePartitionDown
        ( tT,  t0,
               tau1,
         /**/ /****/
          tB,  t2 );

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX
