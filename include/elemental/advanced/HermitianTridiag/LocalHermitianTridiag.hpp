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

namespace {

template<typename R> // representation of a real number
void
HermitianTridiagL( Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("HermitianTridiagL");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
#endif
    // Matrix views 
    Matrix<R>
        ATL, ATR,  A00, a01,     A02,  alpha21T,
        ABL, ABR,  a10, alpha11, a12,  a21B,
                   A20, a21,     A22;

    // Temporary matrices
    Matrix<R> w21;

    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height()+1 < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        PartitionDown
        ( a21, alpha21T,
               a21B,     1 );

        w21.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//
        R tau = advanced::Reflector( alpha21T, a21B );

        R epsilon1 = alpha21T.Get(0,0);
        alpha21T.Set(0,0,(R)1);

        basic::Symv( LOWER, tau, A22, a21, (R)0, w21 );
        R alpha = -static_cast<R>(0.5)*tau*basic::Dot( w21, a21 );
        basic::Axpy( alpha, a21, w21 );
        basic::Syr2( LOWER, (R)-1, a21, w21, A22 );

        alpha21T.Set(0,0,epsilon1);
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

template<typename R> // representation of a real number
void
HermitianTridiagU( Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("HermitianTridiagU");
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
#endif
    // Matrix views 
    Matrix<R>
        ATL, ATR,  A00, a01,     A02,  a01T,
        ABL, ABR,  a10, alpha11, a12,  alpha01B,
                   A20, a21,     A22;

    // Temporary matrices
    Matrix<R> w01;

    PushBlocksizeStack( 1 );
    PartitionUpDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ABR.Height()+1 < A.Height() )
    {
        RepartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );

        PartitionUp
        ( a01, a01T,
               alpha01B, 1 );

        w01.ResizeTo( a01.Height(), 1 );
        //--------------------------------------------------------------------//
        R tau = advanced::Reflector( alpha01B, a01T );

        R epsilon1 = alpha01B.Get(0,0);
        alpha01B.Set(0,0,(R)1);

        basic::Symv( UPPER, tau, A00, a01, (R)0, w01 );
        R alpha = -static_cast<R>(0.5)*tau*basic::Dot( w01, a01 );
        basic::Axpy( alpha, a01, w01 );
        basic::Syr2( UPPER, (R)-1, a01, w01, A00 );

        alpha01B.Set(0,0,epsilon1);
        //--------------------------------------------------------------------//

        SlidePartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
void
HermitianTridiagL
( Matrix< complex<R> >& A, Matrix< complex<R> >& t )
{
#ifndef RELEASE
    PushCallStack("HermitianTridiagL");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( t.Viewing() && (t.Height() != A.Height()-1 || t.Width() != 1) )
        throw std::logic_error
        ("t must be a vector of the same height as A minus one");
#endif
    typedef complex<R> C;

    if( !t.Viewing() )
        t.ResizeTo( A.Height()-1, 1 );

    // Matrix views 
    Matrix<C>
        ATL, ATR,  A00, a01,     A02,  alpha21T,
        ABL, ABR,  a10, alpha11, a12,  a21B,
                   A20, a21,     A22;
    Matrix<C>
        tT,  t0,
        tB,  tau1,
             t2;

    // Temporary matrices
    Matrix<C> w21;

    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( t, tT,
         tB, 0 );
    while( ATL.Height()+1 < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        RepartitionDown
        ( tT,  t0, 
         /**/ /**/
               tau1,
          tB,  t2 );

        PartitionDown
        ( a21, alpha21T,
               a21B,     1 );

        w21.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//
        C tau = advanced::Reflector( alpha21T, a21B );
        tau1.Set(0,0,tau);

        R epsilon1 = alpha21T.GetReal(0,0);
        alpha21T.Set(0,0,(C)1);

        basic::Hemv( LOWER, tau, A22, a21, (C)0, w21 );
        C alpha = -static_cast<C>(0.5)*tau*basic::Dot( w21, a21 );
        basic::Axpy( alpha, a21, w21 );
        basic::Her2( LOWER, (C)-1, a21, w21, A22 );

        alpha21T.Set(0,0,epsilon1);
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );

        SlidePartitionDown
        ( tT,  t0,
               tau1,
         /**/ /****/
          tB,  t2 );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> // representation of a real number
void
HermitianTridiagU
( Matrix< complex<R> >& A, Matrix< complex<R> >& t )
{
#ifndef RELEASE
    PushCallStack("HermitianTridiagU");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( t.Viewing() && (t.Height() != A.Height()-1 || t.Width() != 1) )
        throw std::logic_error
        ("t must be a vector of the same height as A minus one");
#endif
    typedef complex<R> C;
    
    if( !t.Viewing() )
        t.ResizeTo( A.Height()-1, 1 );

    // Matrix views 
    Matrix<C>
        ATL, ATR,  A00, a01,     A02,  a01T,
        ABL, ABR,  a10, alpha11, a12,  alpha01B,
                   A20, a21,     A22;
    Matrix<C>
        tT,  t0,
        tB,  tau1,
             t2;

    // Temporary matrices
    Matrix<C> w01;

    PushBlocksizeStack( 1 );
    PartitionUpDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionUp
    ( t, tT,
         tB, 0 );
    while( ABR.Height()+1 < A.Height() )
    {
        RepartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );

        RepartitionUp
        ( tT,  t0,
               tau1,
         /**/ /****/
          tB,  t2 );

        PartitionUp
        ( a01, a01T,
               alpha01B, 1 );

        w01.ResizeTo( a01.Height(), 1 );
        //--------------------------------------------------------------------//
        C tau = advanced::Reflector( alpha01B, a01T );
        tau1.Set(0,0,tau);

        R epsilon1 = alpha01B.GetReal(0,0);
        alpha01B.Set(0,0,(C)1);

        basic::Hemv( UPPER, tau, A00, a01, (C)0, w01 );
        C alpha = -static_cast<C>(0.5)*tau*basic::Dot( w01, a01 );
        basic::Axpy( alpha, a01, w01 );
        basic::Her2( UPPER, (C)-1, a01, w01, A00 );

        alpha01B.Set(0,0,epsilon1);
        //--------------------------------------------------------------------//

        SlidePartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        SlidePartitionUp
        ( tT,  t0,
         /**/ /****/
               tau1,
          tB,  t2 );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

} // anonymous namespace

template<typename R> // representation of a real number
inline void
elemental::advanced::HermitianTridiag
( Shape shape, Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("HermitianTridiag");
#endif
    if( shape == LOWER )
        HermitianTridiagL( A );
    else
        HermitianTridiagU( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
inline void
elemental::advanced::HermitianTridiag
( Shape shape, Matrix< complex<R> >& A, Matrix< complex<R> >& t )
{
#ifndef RELEASE
    PushCallStack("HermitianTridiag");
#endif
    if( shape == LOWER )
        HermitianTridiagL( A, t );
    else
        HermitianTridiagU( A, t );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX
