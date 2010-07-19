/*
   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Elemental.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "elemental/blas_internal.hpp"
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;

namespace {

template<typename R>
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

template<typename R>
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
template<typename R>
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

template<typename R>
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

template<typename R>
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
template<typename R>
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

