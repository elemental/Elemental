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

namespace elemental {

namespace bidiag {

template<typename R>
inline void BidiagL( Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("bidiag::BidiagL");
    if( A.Height() > A.Width() )
        throw std::logic_error("A must be at least as wide as it is tall");
#endif
    // Matrix views 
    Matrix<R>
        ATL, ATR,  A00, a01,     A02,  alpha21T,  a1R,
        ABL, ABR,  a10, alpha11, a12,  a21B,      A2R,
                   A20, a21,     A22;

    // Temporary matrices
    Matrix<R> w12, w21;

    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        a1R.View1x2( alpha11, a12 );
        A2R.View1x2( a21, A22 );

        w12.ResizeTo( 1, a12.Width() );
        w21.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//
        const R tauP = Reflector( alpha11, a12 );
        const R epsilonP = alpha11.Get(0,0);
        alpha11.Set(0,0,(R)1);
        Gemv( NORMAL, tauP, A2R, a1R, (R)0, w21 );
        Ger( (R)-1, w21, a1R, A2R );
        alpha11.Set(0,0,epsilonP);

        if( A22.Height() != 0 )
        {
            PartitionDown
            ( a21, alpha21T,
                   a21B );

            const R tauQ = Reflector( alpha21T, a21B );
            const R epsilonQ = alpha21T.Get(0,0);
            alpha21T.Set(0,0,(R)1);
            Gemv( TRANSPOSE, tauQ, A22, a21, (R)0, w12 );
            Ger( (R)-1, a21, w12, A22 );
            alpha21T.Set(0,0,epsilonQ);
        }
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
inline void BidiagU( Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("bidiag::BidiagU");
    if( A.Height() < A.Width() )
        throw std::logic_error("A must be at least as tall as it is wide");
#endif
    // Matrix views 
    Matrix<R>
        ATL, ATR,  A00, a01,     A02,  alpha12L, a12R,
        ABL, ABR,  a10, alpha11, a12,  aB1, AB2,
                   A20, a21,     A22;

    // Temporary matrices
    Matrix<R> w12, w21;

    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        aB1.View2x1
        ( alpha11,
          a21 );
        AB2.View2x1
        ( a12,
          A22 );

        w12.ResizeTo( 1, a12.Width() );
        w21.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//
        const R tauQ = Reflector( alpha11, a21 );
        const R epsilonQ = alpha11.Get(0,0);
        alpha11.Set(0,0,(R)1);
        Gemv( TRANSPOSE, tauQ, AB2, aB1, (R)0, w12 );
        Ger( (R)-1, aB1, w12, AB2 );
        alpha11.Set(0,0,epsilonQ);

        if( A22.Width() != 0 )
        {
            PartitionRight( a12, alpha12L, a12R );

            const R tauP = Reflector( alpha12L, a12R );
            const R epsilonP = alpha12L.Get(0,0);
            alpha12L.Set(0,0,(R)1);
            Gemv( NORMAL, tauP, A22, a12, (R)0, w21 );
            Ger( (R)-1, w21, a12, A22 );
            alpha12L.Set(0,0,epsilonP);
        }
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
inline void BidiagL
( Matrix<std::complex<R> >& A,
  Matrix<std::complex<R> >& tP,
  Matrix<std::complex<R> >& tQ )
{
#ifndef RELEASE
    PushCallStack("bidiag::BidiagL");
#endif
    const int tPHeight = A.Height();
    const int tQHeight = std::max(A.Height()-1,0);
#ifndef RELEASE
    if( A.Height() > A.Width() )
        throw std::logic_error("A must be at least as wide as it is tall");
    if( tP.Viewing() && (tP.Height() != tPHeight || tP.Width() != 1) )
        throw std::logic_error("tP is the wrong height");
    if( tQ.Viewing() && (tQ.Height() != tQHeight || tQ.Width() != 1) )
        throw std::logic_error("tQ is the wrong height");
#endif
    typedef std::complex<R> C;

    if( !tP.Viewing() )
        tP.ResizeTo( tPHeight, 1 );
    if( !tQ.Viewing() )
        tQ.ResizeTo( tQHeight, 1 );

    // Matrix views 
    Matrix<C>
        ATL, ATR,  A00, a01,     A02,  alpha21T,  a1R,
        ABL, ABR,  a10, alpha11, a12,  a21B,      A2R,
                   A20, a21,     A22;

    // Temporary matrices
    Matrix<C> w12, w21;

    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        a1R.View1x2( alpha11, a12 );
        A2R.View1x2( a21, A22 );

        w12.ResizeTo( 1, a12.Width() );
        w21.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//
        Conjugate( a1R );
        const C tauP = Reflector( alpha11, a12 );
        const C epsilonP = alpha11.Get(0,0);
        alpha11.Set(0,0,(C)1);
        Gemv( NORMAL, tauP, A2R, a1R, (C)0, w21 );
        Ger( (C)-1, w21, a1R, A2R );
        alpha11.Set(0,0,epsilonP);
        tP.Set(A00.Height(),0,tauP);
        Conjugate( a1R );

        if( A22.Height() != 0 )
        {
            PartitionDown
            ( a21, alpha21T,
                   a21B );

            const C tauQ = Reflector( alpha21T, a21B );
            const C epsilonQ = alpha21T.Get(0,0);
            alpha21T.Set(0,0,(C)1);
            Gemv( ADJOINT, tauQ, A22, a21, (C)0, w12 );
            Ger( (C)-1, a21, w12, A22 );
            alpha21T.Set(0,0,epsilonQ);
            tQ.Set(A00.Height(),0,tauQ);
        }
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
inline void BidiagU
( Matrix<std::complex<R> >& A, 
  Matrix<std::complex<R> >& tP,
  Matrix<std::complex<R> >& tQ )
{
#ifndef RELEASE
    PushCallStack("BidiagU");
#endif
    const int tPHeight = std::max(A.Width()-1,0);
    const int tQHeight = A.Width();
#ifndef RELEASE
    if( A.Height() < A.Width() )
        throw std::logic_error("A must be at least as tall as it is wide");
    if( tP.Viewing() && (tP.Height() != tPHeight || tP.Width() != 1) )
        throw std::logic_error("tP is the wrong height");
    if( tQ.Viewing() && (tQ.Height() != tQHeight || tQ.Width() != 1) )
        throw std::logic_error("tQ is the wrong height");
#endif
    typedef std::complex<R> C;

    if( !tP.Viewing() )
        tP.ResizeTo( tPHeight, 1 );
    if( !tQ.Viewing() )
        tQ.ResizeTo( tQHeight, 1 );

    // Matrix views 
    Matrix<C>
        ATL, ATR,  A00, a01,     A02,  alpha12L, a12R,
        ABL, ABR,  a10, alpha11, a12,  aB1, AB2,
                   A20, a21,     A22;

    // Temporary matrices
    Matrix<C> w12, w21;

    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        aB1.View2x1
        ( alpha11,
          a21 );
        AB2.View2x1
        ( a12,
          A22 );

        w12.ResizeTo( 1, a12.Width() );
        w21.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//
        const C tauQ = Reflector( alpha11, a21 );
        const C epsilonQ = alpha11.Get(0,0);
        alpha11.Set(0,0,(C)1);
        Gemv( ADJOINT, tauQ, AB2, aB1, (C)0, w12 );
        Ger( (C)-1, aB1, w12, AB2 );
        alpha11.Set(0,0,epsilonQ);
        tQ.Set(A00.Height(),0,tauQ );

        if( A22.Width() != 0 )
        {
            PartitionRight( a12, alpha12L, a12R );

            Conjugate( a12 ); 
            const C tauP = Reflector( alpha12L, a12R );
            const C epsilonP = alpha12L.Get(0,0);
            alpha12L.Set(0,0,(C)1);
            Gemv( NORMAL, tauP, A22, a12, (C)0, w21 );
            Ger( (C)-1, w21, a12, A22 );
            alpha12L.Set(0,0,epsilonP);
            Conjugate( a12 );
            tP.Set(A00.Height(),0,tauP);
        }
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

} // namespace bidiag

template<typename R> 
inline void
Bidiag( Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("Bidiag");
#endif
    if( IsComplex<R>::val )
        throw std::logic_error("Called real routine with complex datatype");
    if( A.Height() >= A.Width() )
        bidiag::BidiagU( A );
    else
        bidiag::BidiagL( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void Bidiag
( Matrix<std::complex<R> >& A, 
  Matrix<std::complex<R> >& tP, 
  Matrix<std::complex<R> >& tQ )
{
#ifndef RELEASE
    PushCallStack("Bidiag");
#endif
    if( A.Height() >= A.Width() )
        bidiag::BidiagU( A, tP, tQ );
    else
        bidiag::BidiagL( A, tP, tQ );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elemental
