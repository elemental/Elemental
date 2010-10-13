/*
   Copyright (c) 2009-2010, Jack Poulson
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

template<typename R>
void
HessenbergL
( Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("HessenbergL");
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
#endif
    // Matrix views 
    Matrix<R>
        ATL, ATR,  A00, a01,     A02,  alpha21T,
        ABL, ABR,  a10, alpha11, a12,  a21B,
                   A20, a21,     A22;

    // Temporary matrices
    Matrix<R> w01;
    Matrix<R> y21;
    Matrix<R> z21;

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

        PartitionDown
        ( a21, alpha21T,
               a21B,     1 );

        w01.ResizeTo( a01.Height(), 1 );
        y21.ResizeTo( a21.Height(), 1 );
        z21.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//
        R tau = lapack::Reflector( alpha21T, a21B );
        R epsilon = alpha21T(0,0);
        alpha21T(0,0) = (R)1;

        blas::Gemv( Transpose, (R)1, A22, a21, (R)0, y21 );
        blas::Gemv( Normal,    (R)1, A22, a21, (R)0, z21 );

        R beta = 0.5*blas::Dot( a21, z21 );
        blas::Scal( tau, y21 );
        blas::Scal( tau, z21 );
        blas::Axpy( -beta*tau, a21, y21 );
        blas::Axpy( -beta*tau, a21, z21 );

        beta = blas::Dot( a12, a21 );
        blas::Axpy( -beta*tau, a21, a12 );

        blas::Gemv( Normal, (R)1, A02, a21, (R)0, w01 );
        blas::Ger( -tau,  w01, a21, A02 );
        blas::Ger( (R)-1, a21, y21, A22 );
        blas::Ger( (R)-1, z21, a21, A22 );

        alpha21T(0,0) = epsilon;
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
void
HessenbergU
( Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("HessenbergU");
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
#endif
    // Matrix views 
    Matrix<R>
        ATL, ATR,  A00, a01,     A02,  a01T,
        ABL, ABR,  a10, alpha11, a12,  alpha01B,
                   A20, a21,     A22;

    // Temporary matrices
    Matrix<R> y01;
    Matrix<R> z01;
    Matrix<R> w21;

    PushBlocksizeStack( 1 );
    PartitionUpDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ABR.Height() < A.Height() )
    {
        RepartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );

        PartitionUp
        ( a01, a01T,
               alpha01B, 1 );

        y01.ResizeTo( a01.Height(), 1 );
        z01.ResizeTo( a01.Height(), 1 );
        w21.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//
        R tau = lapack::Reflector( alpha01B, a01T );
        R epsilon1 = alpha01B(0,0);
        alpha01B(0,0) = (R)1;

        blas::Gemv( Transpose, (R)1, A00, a01, (R)0, y01 );
        blas::Gemv( Normal,    (R)1, A00, a01, (R)0, z01 );

        R beta = 0.5*blas::Dot( a01, z01 );
        blas::Scal( tau, y01 );
        blas::Scal( tau, z01 );
        blas::Axpy( -beta*tau, a01, y01 );
        blas::Axpy( -beta*tau, a01, z01 );

        blas::Dot( a10, a01 );
        blas::Axpy( -beta*tau, a01, a10 );

        blas::Gemv( Normal, (R)1, A20, a01, (R)0, w21 );
        blas::Ger( -tau,  w21, a01, A20 );
        blas::Ger( (R)-1, a01, y01, A00 );
        blas::Ger( (R)-1, z01, a01, A00 );

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
template<typename R>
void
HessenbergL
( Matrix< complex<R> >& A, Matrix< complex<R> >& t )
{
#ifndef RELEASE
    PushCallStack("TridiagL");
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
    if( t.Viewing() && (t.Height() != A.Height() || t.Width() != 1) )
        throw logic_error
              ( "t must be a vector of the same height as A." );
#endif
    typedef complex<R> C;

    if( !t.Viewing() )
        t.ResizeTo( A.Height(), 1 );

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
    Matrix<C> w01;
    Matrix<C> y21;
    Matrix<C> z21;

    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( t, tT,
         tB, 0 );
    while( ATL.Height() < A.Height() )
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

        w01.ResizeTo( a01.Height(), 1 );
        y21.ResizeTo( a21.Height(), 1 );
        z21.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//
        C tau = lapack::Reflector( alpha21T, a21B );
        tau1(0,0) = tau;
        C epsilon1 = alpha21T(0,0);
        alpha21T(0,0) = (C)1;

        blas::Gemv( ConjugateTranspose, (C)1, A22, a21, (C)0, y21 );
        blas::Gemv( Normal,             (C)1, A22, a21, (C)0, z21 );

        C beta = static_cast<C>(0.5)*blas::Dot( a21, z21 );
        blas::Scal( tau, y21 );
        blas::Scal( tau, z21 );
        blas::Axpy( -beta*tau, a21, y21 );
        blas::Axpy( -beta*tau, a21, z21 );

        beta = blas::Dot( a12, a21 );
        blas::Axpy( -beta*tau, a21, a12 );

        blas::Gemv( Normal, (C)1, A02, a21, (C)0, w01 );
        blas::Ger( -tau, w01, a21, A02 );
        blas::Ger( (C)-1, a21, y21, A22 );
        blas::Ger( (C)-1, z21, a21, A22 );

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

template<typename R>
void
HessenbergU
( Matrix< complex<R> >& A, Matrix< complex<R> >& t )
{
#ifndef RELEASE
    PushCallStack("TridiagU");
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
    if( t.Viewing() && (t.Height() != A.Height() || t.Width() != 1) )
        throw logic_error
              ( "t must be a vector of the same height as A." );
#endif
    typedef complex<R> C;
    
    if( !t.Viewing() )
        t.ResizeTo( A.Height(), 1 );

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
    Matrix<C> y01;
    Matrix<C> z01;
    Matrix<C> w21;

    PushBlocksizeStack( 1 );
    PartitionUpDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionUp
    ( t, tT,
         tB, 0 );
    while( ABR.Height() < A.Height() )
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

        y01.ResizeTo( a01.Height(), 1 );
        z01.ResizeTo( a01.Height(), 1 );
        w21.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//
        C tau = lapack::Reflector( alpha01B, a01T );
        tau1(0,0) = tau;
        C epsilon1 = alpha01B(0,0);
        alpha01B(0,0) = (C)1;

        blas::Gemv( ConjugateTranspose, (C)1, A00, a01, (C)0, y01 );
        blas::Gemv( Normal,             (C)1, A00, a01, (C)0, z01 );

        C beta = static_cast<C>(0.5)*blas::Dot( a01, z01 );
        blas::Scal( tau, y01 );
        blas::Scal( tau, z01 );
        blas::Axpy( -beta*tau, a01, y01 );
        blas::Axpy( -beta*tau, a01, z01 );

        blas::Dot( a10, a01 );
        blas::Axpy( -beta*tau, a01, a10 );

        blas::Gemv( Normal, (C)1, A20, a01, (C)0, w21 );
        blas::Ger( -tau,  w21, a01, A20 );
        blas::Ger( (C)-1, a01, y01, A00 );
        blas::Ger( (C)-1, z01, a01, A00 );

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

template<typename R>
void
elemental::lapack::Hessenberg
( Shape shape, Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("Hessenberg");
#endif
    if( shape == Lower )
        HessenbergL( A );
    else
        HessenbergU( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::lapack::Hessenberg
( Shape shape, Matrix< complex<R> >& A, Matrix< complex<R> >& t )
{
#ifndef RELEASE
    PushCallStack("Hessenberg");
#endif
    if( shape == Lower )
        HessenbergL( A, t );
    else
        HessenbergU( A, t );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template void elemental::lapack::Hessenberg
( Shape shape, Matrix<float>& A );

template void elemental::lapack::Hessenberg
( Shape shape, Matrix<double>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::Hessenberg
( Shape shape, Matrix<scomplex>& A, Matrix<scomplex>& t );

template void elemental::lapack::Hessenberg
( Shape shape, Matrix<dcomplex>& A, Matrix<dcomplex>& t );
#endif

