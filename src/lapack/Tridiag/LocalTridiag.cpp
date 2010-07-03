/*
   This file is part of Elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (c) 2009-2010 Jack Poulson <jack.poulson@gmail.com>.
   All rights reserved.

   This file is released under the terms of the license contained in the file
   LICENSE-PURE.
*/
#include "elemental/blas_internal.hpp"
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;

namespace {

template<typename R>
R
Reflector
( Matrix<R>& chi, Matrix<R>& x )
{
#ifndef RELEASE
    PushCallStack("Reflector");
#endif
    if( x.Height() == 0 )
    {
        chi(0,0) *= (R)-1;
#ifndef RELEASE
        PopCallStack();
#endif
        return (R)2;
    }

    R norm = blas::Nrm2( x );
    R alpha = chi(0,0);
    
    R beta;
    if( alpha <= 0 )
        beta = wrappers::lapack::SafeNorm( alpha, norm );
    else
        beta = -wrappers::lapack::SafeNorm( alpha, norm );

    R safeMin = numeric_limits<R>::min() / numeric_limits<R>::epsilon();
    int count = 0;
    if( Abs( beta ) < safeMin )
    {
        R invOfSafeMin = static_cast<R>(1) / safeMin;
        do
        {
            ++count;
            blas::Scal( invOfSafeMin, x );
            alpha *= invOfSafeMin;
            beta *= invOfSafeMin;
        } while( Abs( beta ) < safeMin );

        norm = blas::Nrm2( x );
        if( alpha <= 0 )
            beta = wrappers::lapack::SafeNorm( alpha, norm );
        else
            beta = -wrappers::lapack::SafeNorm( alpha, norm );
    }

    R tau = ( beta - alpha ) / beta;
    blas::Scal( static_cast<R>(1)/(alpha-beta), x );

    for( int j=0; j<count; ++j )
        beta *= safeMin;
    chi(0,0) = beta;
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}

#ifndef WITHOUT_COMPLEX
template<typename R>
complex<R>
Reflector
( Matrix< complex<R> >& chi, Matrix< complex<R> >& x )
{
#ifndef RELEASE
    PushCallStack("Reflector");
#endif
    typedef complex<R> C;

    R norm = blas::Nrm2( x );
    C alpha = chi(0,0);

    if( norm == 0 && imag(alpha) == (R)0 )
    {
        chi(0,0) *= (R)-1;
#ifndef RELEASE
        PopCallStack();
#endif
        return (C)2;
    }

    R beta;
    if( real(alpha) <= 0 )
        beta = wrappers::lapack::SafeNorm( real(alpha), imag(alpha), norm );
    else
        beta = -wrappers::lapack::SafeNorm( real(alpha), imag(alpha), norm );

    R safeMin = numeric_limits<R>::min() / numeric_limits<R>::epsilon();
    int count = 0;
    if( Abs( beta ) < safeMin )
    {
        R invOfSafeMin = static_cast<R>(1) / safeMin;
        do
        {
            ++count;
            blas::Scal( (C)invOfSafeMin, x );
            alpha *= invOfSafeMin;
            beta *= invOfSafeMin;
        } while( Abs( beta ) < safeMin );

        norm = blas::Nrm2( x );
        if( real(alpha) <= 0 )
            beta = wrappers::lapack::SafeNorm
                   ( real(alpha), imag(alpha), norm );
        else
            beta = -wrappers::lapack::SafeNorm
                    ( real(alpha), imag(alpha), norm );
    }

    C tau = C( (beta-real(alpha))/beta, -imag(alpha)/beta );
    blas::Scal( static_cast<C>(1)/(alpha-beta), x );

    for( int j=0; j<count; ++j )
        beta *= safeMin;
    chi(0,0) = beta;
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}
#endif // WITHOUT_COMPLEX

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
    Matrix<R> z;

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

        z.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//
        R tau = Reflector( alpha21T, a21B );

        R epsilon1 = alpha21T(0,0);
        alpha21T(0,0) = (R)1;

        blas::Symv( Lower, (R)1, A22, a21, (R)0, z );
        R alpha = -static_cast<R>(0.5)*tau*blas::Dot( z, a21 );
        blas::Axpy( alpha, a21, z );
        blas::Syr2( Lower, (R)-1, a21, z, A22 );

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
    Matrix<R> z;

    PushBlocksizeStack( 1 );
    PartitionUpDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ABR.Height()+1 < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );

        PartitionUp
        ( a01, a01T,
               alpha01B, 1 );

        z.ResizeTo( a01.Height(), 1 );
        //--------------------------------------------------------------------//
        R tau = Reflector( alpha01B, a01T );

        R epsilon1 = alpha01B(0,0);
        alpha01B(0,0) = (R)1;

        blas::Symv( Upper, (R)1, A22, a01, (R)0, z );
        R alpha = -static_cast<R>(0.5)*tau*blas::Dot( z, a01 );
        blas::Axpy( alpha, a01, z );
        blas::Syr2( Upper, (R)-1, a01, z, A22 );

        alpha01B(0,0) = epsilon1;
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
    if( t.Height() != A.Height()-1 || t.Width() != 1 )
        throw logic_error
              ( "t must be a vector of the same height as A minus one." );
#endif
    typedef complex<R> C;

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
    Matrix<C> z;

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

        z.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//
        C tau = Reflector( alpha21T, a21B );
        tau1(0,0) = tau;

        R epsilon1 = real(alpha21T(0,0));
        alpha21T(0,0) = (C)1;

        blas::Hemv( Lower, (C)1, A22, a21, (C)0, z );
        C alpha = -static_cast<C>(0.5)*tau*blas::Dot( z, a21 );
        blas::Axpy( alpha, a21, z );
        blas::Her2( Lower, (C)-1, a21, z, A22 );

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
    if( t.Height() != A.Height()-1 || t.Width() != 1 )
        throw logic_error
              ( "t must be a vector of the same height as A minus one." );
#endif
    typedef complex<R> C;

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
    Matrix<C> z;

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

        z.ResizeTo( a01.Height(), 1 );
        //--------------------------------------------------------------------//
        C tau = Reflector( alpha01B, a01T );
        tau1(0,0) = tau;

        R epsilon1 = real(alpha01B(0,0));
        alpha01B(0,0) = (C)1;

        blas::Hemv( Upper, (C)1, A22, a01, (C)0, z );
        C alpha = -static_cast<C>(0.5)*tau*blas::Dot( z, a01 );
        blas::Axpy( alpha, a01, z );
        blas::Her2( Upper, (C)-1, a01, z, A22 );

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

