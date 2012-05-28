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

template<typename R>
inline void
internal::PanelQR( Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("internal::PanelQR");
#endif
    Matrix<R>
        ATL, ATR,  A00, a01,     A02,  aLeftCol, ARightPan,
        ABL, ABR,  a10, alpha11, a12,
                   A20, a21,     A22;

    Matrix<R> z;

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

        Zeros( ARightPan.Width(), 1, z );
        //--------------------------------------------------------------------//
        const R tau = Reflector( alpha11, a21 );
        const R alpha = alpha11.Get(0,0);
        alpha11.Set(0,0,1);

        Gemv( TRANSPOSE, (R)1, ARightPan, aLeftCol, (R)0, z );
        Ger( -tau, aLeftCol, z, ARightPan );

        alpha11.Set(0,0,alpha);
        //--------------------------------------------------------------------//

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

template<typename R>
inline void
internal::PanelQR( DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("internal::PanelQR");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<R,MC,MR>
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  aLeftCol(g), ARightPan(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),
                         A20(g), a21(g),     A22(g);

    // Temporary distributions
    DistMatrix<R,MC,STAR> aLeftCol_MC_STAR(g);
    DistMatrix<R,MR,STAR> z_MR_STAR(g);

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
        z_MR_STAR.AlignWith( ARightPan );
        Zeros( ARightPan.Width(), 1, z_MR_STAR );
        //--------------------------------------------------------------------//
        const R tau = Reflector( alpha11, a21 );

        const bool myDiagonalEntry = ( g.Row() == alpha11.ColAlignment() && 
                                       g.Col() == alpha11.RowAlignment() );
        R alpha = (R)0;
        if( myDiagonalEntry )
        {
            alpha = alpha11.GetLocalEntry(0,0);
            alpha11.SetLocalEntry(0,0,1);
        }

        aLeftCol_MC_STAR = aLeftCol;

        Gemv
        ( TRANSPOSE, 
          (R)1, ARightPan.LockedLocalMatrix(), 
                aLeftCol_MC_STAR.LockedLocalMatrix(),
          (R)0, z_MR_STAR.LocalMatrix() );
        z_MR_STAR.SumOverCol(); 

        Ger
        ( -tau, 
          aLeftCol_MC_STAR.LockedLocalMatrix(), 
          z_MR_STAR.LockedLocalMatrix(),
          ARightPan.LocalMatrix() );

        if( myDiagonalEntry )
            alpha11.SetLocalEntry(0,0,alpha);
        //--------------------------------------------------------------------//
        aLeftCol_MC_STAR.FreeAlignments();
        z_MR_STAR.FreeAlignments();

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

template<typename R> 
inline void
internal::PanelQR
( Matrix<Complex<R> >& A,
  Matrix<Complex<R> >& t )
{
#ifndef RELEASE
    PushCallStack("internal::PanelQR");
    if( t.Height() != std::min(A.Height(),A.Width()) || t.Width() != 1 )
        throw std::logic_error
        ("t must be a vector of height equal to the minimum dimension of A");
#endif
    typedef Complex<R> C;

    Matrix<C>
        ATL, ATR,  A00, a01,     A02,  aLeftCol, ARightPan,
        ABL, ABR,  a10, alpha11, a12,
                   A20, a21,     A22;
    DistMatrix<C,MD,STAR>
        tT,  t0,
        tB,  tau1,
             t2;

    Matrix<C> z;

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

        Zeros( ARightPan.Width(), 1, z );
        //--------------------------------------------------------------------//
        const C tau = Reflector( alpha11, a21 );
        tau1.Set( 0, 0, tau );
        const C alpha = alpha11.Get(0,0);
        alpha11.Set(0,0,1);

        Gemv( ADJOINT, (C)1, ARightPan, aLeftCol, (C)0, z );
        Ger( -Conj(tau), aLeftCol, z, ARightPan );

        alpha11.Set(0,0,alpha);
        //--------------------------------------------------------------------//

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

template<typename R> 
inline void
internal::PanelQR
( DistMatrix<Complex<R>,MC,MR  >& A,
  DistMatrix<Complex<R>,MD,STAR>& t )
{
#ifndef RELEASE
    PushCallStack("internal::PanelQR");
    if( A.Grid() != t.Grid() )
        throw std::logic_error("{A,t} must be distributed over the same grid");
    if( t.Height() != std::min(A.Height(),A.Width()) || t.Width() != 1 )
        throw std::logic_error
        ("t must be a vector of height equal to the minimum dimension of A");
    if( !t.AlignedWithDiagonal( A, 0 ) )
        throw std::logic_error("t must be aligned with A's main diagonal");
#endif
    typedef Complex<R> C;
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
    DistMatrix<C,MR,STAR> z_MR_STAR(g);

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
        z_MR_STAR.AlignWith( ARightPan );
        Zeros( ARightPan.Width(), 1, z_MR_STAR );
        //--------------------------------------------------------------------//
        const C tau = Reflector( alpha11, a21 );
        tau1.Set( 0, 0, tau );

        const bool myDiagonalEntry = ( g.Row() == alpha11.ColAlignment() && 
                                       g.Col() == alpha11.RowAlignment() );
        C alpha = (C)0;
        if( myDiagonalEntry )
        {
            alpha = alpha11.GetLocalEntry(0,0);
            alpha11.SetLocalEntry(0,0,1);
        }

        aLeftCol_MC_STAR = aLeftCol;

        Gemv
        ( ADJOINT, 
          (C)1, ARightPan.LockedLocalMatrix(), 
                aLeftCol_MC_STAR.LockedLocalMatrix(),
          (C)0, z_MR_STAR.LocalMatrix() );
        z_MR_STAR.SumOverCol(); 

        Ger
        ( -Conj(tau), 
          aLeftCol_MC_STAR.LockedLocalMatrix(), 
          z_MR_STAR.LockedLocalMatrix(),
          ARightPan.LocalMatrix() );

        if( myDiagonalEntry )
            alpha11.SetLocalEntry(0,0,alpha);
        //--------------------------------------------------------------------//
        aLeftCol_MC_STAR.FreeAlignments();
        z_MR_STAR.FreeAlignments();

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

} // namespace elem
