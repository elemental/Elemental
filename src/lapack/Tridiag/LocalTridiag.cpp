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
#include "elemental/blas_internal.hpp"
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;

namespace {

template<typename R> // representation of a real number
void
TridiagL
( Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("TridiagL");
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
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
        R tau = lapack::Reflector( alpha21T, a21B );

        R epsilon1 = alpha21T(0,0);
        alpha21T(0,0) = (R)1;

        blas::Symv( Lower, tau, A22, a21, (R)0, w21 );
        R alpha = -static_cast<R>(0.5)*tau*blas::Dot( w21, a21 );
        blas::Axpy( alpha, a21, w21 );
        blas::Syr2( Lower, (R)-1, a21, w21, A22 );

        alpha21T(0,0) = epsilon1;
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
TridiagU
( Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("TridiagU");
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
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
        R tau = lapack::Reflector( alpha01B, a01T );

        R epsilon1 = alpha01B(0,0);
        alpha01B(0,0) = (R)1;

        blas::Symv( Upper, tau, A00, a01, (R)0, w01 );
        R alpha = -static_cast<R>(0.5)*tau*blas::Dot( w01, a01 );
        blas::Axpy( alpha, a01, w01 );
        blas::Syr2( Upper, (R)-1, a01, w01, A00 );

        alpha01B(0,0) = epsilon1;
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
TridiagL
( Matrix< complex<R> >& A, Matrix< complex<R> >& t )
{
#ifndef RELEASE
    PushCallStack("TridiagL");
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
    if( t.Viewing() && (t.Height() != A.Height()-1 || t.Width() != 1) )
        throw logic_error
              ( "t must be a vector of the same height as A minus one." );
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
        C tau = lapack::Reflector( alpha21T, a21B );
        tau1(0,0) = tau;

        R epsilon1 = real(alpha21T(0,0));
        alpha21T(0,0) = (C)1;

        blas::Hemv( Lower, tau, A22, a21, (C)0, w21 );
        C alpha = -static_cast<C>(0.5)*tau*blas::Dot( w21, a21 );
        blas::Axpy( alpha, a21, w21 );
        blas::Her2( Lower, (C)-1, a21, w21, A22 );

        alpha21T(0,0) = epsilon1;
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
TridiagU
( Matrix< complex<R> >& A, Matrix< complex<R> >& t )
{
#ifndef RELEASE
    PushCallStack("TridiagU");
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
    if( t.Viewing() && (t.Height() != A.Height()-1 || t.Width() != 1) )
        throw logic_error
              ( "t must be a vector of the same height as A minus one." );
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
        C tau = lapack::Reflector( alpha01B, a01T );
        tau1(0,0) = tau;

        R epsilon1 = real(alpha01B(0,0));
        alpha01B(0,0) = (C)1;

        blas::Hemv( Upper, tau, A00, a01, (C)0, w01 );
        C alpha = -static_cast<C>(0.5)*tau*blas::Dot( w01, a01 );
        blas::Axpy( alpha, a01, w01 );
        blas::Her2( Upper, (C)-1, a01, w01, A00 );

        alpha01B(0,0) = epsilon1;
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
void
elemental::lapack::Tridiag
( Shape shape, Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("Tridiag");
#endif
    if( shape == Lower )
        TridiagL( A );
    else
        TridiagU( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
void
elemental::lapack::Tridiag
( Shape shape, Matrix< complex<R> >& A, Matrix< complex<R> >& t )
{
#ifndef RELEASE
    PushCallStack("Tridiag");
#endif
    if( shape == Lower )
        TridiagL( A, t );
    else
        TridiagU( A, t );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template void elemental::lapack::Tridiag
( Shape shape, Matrix<float>& A );

template void elemental::lapack::Tridiag
( Shape shape, Matrix<double>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::Tridiag
( Shape shape, Matrix<scomplex>& A, Matrix<scomplex>& t );

template void elemental::lapack::Tridiag
( Shape shape, Matrix<dcomplex>& A, Matrix<dcomplex>& t );
#endif

