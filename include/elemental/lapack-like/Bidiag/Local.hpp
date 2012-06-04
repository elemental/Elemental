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
    Matrix<R> x12Trans, w21;

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

        x12Trans.ResizeTo( a12.Width(), 1 );
        w21.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//

        // Find tauP, v, and epsilonP such that
        //     I - tauP | 1 | | 1, v^T | | alpha11 | = | epsilonP |
        //              | v |            |  a12^T  | = |    0     |
        const R tauP = Reflector( alpha11, a12 );
        const R epsilonP = alpha11.Get(0,0);

        // Set a1R^T = | 1 | and form w21 := A2R a1R^T = A2R | 1 |
        //             | v |                                 | v |
        alpha11.Set(0,0,(R)1);
        Gemv( NORMAL, (R)1, A2R, a1R, (R)0, w21 );

        // A2R := A2R - tauP w21 a1R
        //      = A2R - tauP A2R a1R^T a1R
        //      = A2R (I - tauP a1R^T a1R)
        Ger( -tauP, w21, a1R, A2R );

        // Put epsilonP back instead of the temporary value, 1
        alpha11.Set(0,0,epsilonP);

        if( A22.Height() != 0 )
        {
            // Expose the subvector we seek to zero, a21B
            PartitionDown
            ( a21, alpha21T,
                   a21B );

            // Find tauQ, u, and epsilonQ such that
            //     I - tauQ | 1 | | 1, u^T | | alpha21T | = | epsilonQ |
            //              | u |            |   a21B   | = |    0     |
            const R tauQ = Reflector( alpha21T, a21B );
            const R epsilonQ = alpha21T.Get(0,0);

            // Set a21 = | 1 | and form x12^T = (a21^T A22)^T = A22^T a21
            //           | u |
            alpha21T.Set(0,0,(R)1);
            Gemv( TRANSPOSE, (R)1, A22, a21, (R)0, x12Trans );

            // A22 := A22 - tauQ a21 x12
            //      = A22 - tauQ a21 a21^T A22
            //      = (I - tauQ a21 a21^T) A22
            Ger( -tauQ, a21, x12Trans, A22 );

            // Put epsilonQ back instead of the temporary value, 1
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
    Matrix<R> x12Trans, w21;

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

        x12Trans.ResizeTo( a12.Width(), 1 );
        w21.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//

        // Find tauQ, u, and epsilonQ such that
        //     I - tauQ | 1 | | 1, u^T | | alpha11 | = | epsilonQ |
        //              | u |            |   a21   | = |    0     |
        const R tauQ = Reflector( alpha11, a21 );
        const R epsilonQ = alpha11.Get(0,0);

        // Set aB1 = | 1 | and form x12^T := (aB1^T AB2)^T = AB2^T aB1
        //           | u |
        alpha11.Set(0,0,(R)1);
        Gemv( TRANSPOSE, (R)1, AB2, aB1, (R)0, x12Trans );

        // Update AB2 := AB2 - tauQ aB1 x12
        //             = AB2 - tauQ aB1 aB1^T AB2
        //             = (I - tauQ aB1 aB1^T) AB2
        Ger( -tauQ, aB1, x12Trans, AB2 );

        // Put epsilonQ back instead of the temporary value, 1
        alpha11.Set(0,0,epsilonQ);

        if( A22.Width() != 0 )
        {
            // Expose the subvector we seek to zero, a12R
            PartitionRight( a12, alpha12L, a12R );

            // Find tauP, v, and epsilonP such that
            //     I - tauP | 1 | | 1, v^T | | alpha12L | = | epsilonP |
            //              | v |            |  a12R^T  | = |    0     |
            const R tauP = Reflector( alpha12L, a12R );
            const R epsilonP = alpha12L.Get(0,0);

            // Set a12^T = | 1 | and form w21 := A22 a12^T = A22 | 1 |
            //             | v |                                 | v |
            alpha12L.Set(0,0,(R)1);
            Gemv( NORMAL, (R)1, A22, a12, (R)0, w21 );

            // A22 := A22 - tauP w21 a12
            //      = A22 - tauP A22 a12^T a12
            //      = A22 (I - tauP a12^T a12)
            Ger( -tauP, w21, a12, A22 );

            // Put epsilonP back instead of the temporary value, 1
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
( Matrix<Complex<R> >& A,
  Matrix<Complex<R> >& tP,
  Matrix<Complex<R> >& tQ )
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
    typedef Complex<R> C;

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
    Matrix<C> x12Adj, w21;

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

        x12Adj.ResizeTo( a12.Width(), 1 );
        w21.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//

        // Due to deficiencies in the BLAS ?gemv routines, this section is 
        // easier if we temporarily conjugate a1R = | alpha11, a12 |
        Conjugate( a1R );

        // Find tauP, v, and epsilonP such that
        //     I - conj(tauP) | 1 | | 1, v^H | | alpha11 | = | epsilonP | 
        //                    | v |            |   a12^T |   |    0     |
        const C tauP = Reflector( alpha11, a12 );
        const C epsilonP = alpha11.Get(0,0);
        tP.Set(A00.Height(),0,tauP);

        // Set a1R^T = | 1 | and form w21 := A2R a1R^T = A2R | 1 |
        //             | v |                                 | v |
        alpha11.Set(0,0,(C)1);
        Gemv( NORMAL, (C)1, A2R, a1R, (C)0, w21 );

        // A2R := A2R - tauP w21 conj(a1R)
        //      = A2R - tauP A2R a1R^T conj(a1R)
        //      = A22 (I - tauP a1R^T conj(a1R))
        //      = A22 conj(I - conj(tauP) a1R^H a1R)
        // which compensates for the fact that the reflector was generated
        // on the conjugated a1R.
        Ger( -tauP, w21, a1R, A2R );

        // Put epsilonP back instead of the temporary value, 1
        alpha11.Set(0,0,epsilonP);

        // Undo the temporary conjugation
        Conjugate( a1R );

        if( A22.Height() != 0 )
        {
            // Expose the subvector we seek to zero, a21B
            PartitionDown
            ( a21, alpha21T,
                   a21B );

            // Find tauQ, u, and epsilonQ such that
            //     I - conj(tauQ) | 1 | | 1, u^H | | alpha21T | = | epsilonQ |
            //                    | u |            |   a21B   | = |    0     |
            const C tauQ = Reflector( alpha21T, a21B );
            const C epsilonQ = alpha21T.Get(0,0);
            tQ.Set(A00.Height(),0,tauQ);

            // Set a21 = | 1 | and form x12^H = (a21^H A22)^H = A22^H a21
            //           | u |
            alpha21T.Set(0,0,(C)1);
            Gemv( ADJOINT, (C)1, A22, a21, (C)0, x12Adj );

            // A22 := A22 - conj(tauQ) a21 x12 
            //      = A22 - conj(tauQ) a21 a21^H A22
            //      = (I - conj(tauQ) a21 a21^H) A22
            Ger( -Conj(tauQ), a21, x12Adj, A22 );

            // Put epsilonQ back instead of the temporary value, 1
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
inline void BidiagU
( Matrix<Complex<R> >& A, 
  Matrix<Complex<R> >& tP,
  Matrix<Complex<R> >& tQ )
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
    typedef Complex<R> C;

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
    Matrix<C> x12Adj, w21;

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

        x12Adj.ResizeTo( a12.Width(),  1 );
        w21.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//

        // Find tauQ, u, and epsilonQ such that
        //     I - conj(tauQ) | 1 | | 1, u^H | | alpha11 | = | epsilonQ |
        //                    | u |            |    a21  |   |    0     |
        const C tauQ = Reflector( alpha11, a21 );
        const C epsilonQ = alpha11.Get(0,0);
        tQ.Set(A00.Height(),0,tauQ );

        // Set aB1 = | 1 | and form x12^H := (aB1^H AB2)^H = AB2^H aB1
        //           | u |
        alpha11.Set(0,0,(C)1);
        Gemv( ADJOINT, (C)1, AB2, aB1, (C)0, x12Adj );

        // Update AB2 := AB2 - conj(tauQ) aB1 x12
        //             = AB2 - conj(tauQ) aB1 aB1^H AB2 
        //             = (I - conj(tauQ) aB1 aB1^H) AB2
        Ger( -Conj(tauQ), aB1, x12Adj, AB2 );

        // Put epsilonQ back instead of the temporary value, 1
        alpha11.Set(0,0,epsilonQ);

        if( A22.Width() != 0 )
        {
            // Due to the deficiencies in the BLAS ?gemv routines, this section
            // is easier if we temporarily conjugate a12
            Conjugate( a12 ); 

            // Expose the subvector we seek to zero, a12R
            PartitionRight( a12, alpha12L, a12R );

            // Find tauP, v, and epsilonP such that
            //     I - conj(tauP) | 1 | | 1, v^H | | alpha12L | = | epsilonP |
            //                    | v |            |  a12R^T  |   |    0     |
            const C tauP = Reflector( alpha12L, a12R );
            const C epsilonP = alpha12L.Get(0,0);
            tP.Set(A00.Height(),0,tauP);

            // Set a12^T = | 1 | and form w21 := A22 a12^T = A22 | 1 |
            //             | v |                                 | v |
            alpha12L.Set(0,0,(C)1);
            Gemv( NORMAL, (C)1, A22, a12, (C)0, w21 );

            // A22 := A22 - tauP w21 conj(a12)
            //      = A22 - tauP A22 a12^T conj(a12)
            //      = A22 (I - tauP a12^T conj(a12))
            //      = A22 conj(I - conj(tauP) a12^H a12)
            // which compensates for the fact that the reflector was generated
            // on the conjugated a12.
            Ger( -tauP, w21, a12, A22 );

            // Put epsilonP back instead of the temporary value, 1
            alpha12L.Set(0,0,epsilonP);

            // Undue the temporary conjugation
            Conjugate( a12 );
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
( Matrix<Complex<R> >& A, 
  Matrix<Complex<R> >& tP, 
  Matrix<Complex<R> >& tQ )
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

} // namespace elem
