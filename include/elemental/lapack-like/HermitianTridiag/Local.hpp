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

namespace hermitian_tridiag {

template<typename R> 
inline void HermitianTridiagL( Matrix<R>& A )
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
        const R tau = Reflector( alpha21T, a21B );
        const R epsilon1 = alpha21T.Get(0,0);
        alpha21T.Set(0,0,R(1));

        Symv( LOWER, tau, A22, a21, R(0), w21 );
        const R alpha = -tau*Dot( w21, a21 )/R(2);
        Axpy( alpha, a21, w21 );
        Syr2( LOWER, R(-1), a21, w21, A22 );
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

template<typename R> 
inline void HermitianTridiagU( Matrix<R>& A )
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
        const R tau = Reflector( alpha01B, a01T );
        const R epsilon1 = alpha01B.Get(0,0);
        alpha01B.Set(0,0,R(1));

        Symv( UPPER, tau, A00, a01, R(0), w01 );
        const R alpha = -tau*Dot( w01, a01 )/R(2);
        Axpy( alpha, a01, w01 );
        Syr2( UPPER, R(-1), a01, w01, A00 );
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

template<typename R>
inline void HermitianTridiagL
( Matrix<Complex<R> >& A, Matrix<Complex<R> >& t )
{
#ifndef RELEASE
    PushCallStack("HermitianTridiagL");
#endif
    const int tHeight = std::max(A.Height()-1,0);
#ifndef RELEASE
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( t.Viewing() && (t.Height() != tHeight || t.Width() != 1) )
        throw std::logic_error("t is of the wrong size");
#endif
    typedef Complex<R> C;

    if( !t.Viewing() )
        t.ResizeTo( tHeight, 1 );

    // Matrix views 
    Matrix<C>
        ATL, ATR,  A00, a01,     A02,  alpha21T,
        ABL, ABR,  a10, alpha11, a12,  a21B,
                   A20, a21,     A22;

    // Temporary matrices
    Matrix<C> w21;

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
        const C tau = Reflector( alpha21T, a21B );
        const R epsilon1 = alpha21T.GetRealPart(0,0);
        t.Set(A00.Height(),0,tau);
        alpha21T.Set(0,0,C(1));

        Hemv( LOWER, tau, A22, a21, C(0), w21 );
        const C alpha = -tau*Dot( w21, a21 )/C(2);
        Axpy( alpha, a21, w21 );
        Her2( LOWER, C(-1), a21, w21, A22 );
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

template<typename R>
inline void HermitianTridiagU
( Matrix<Complex<R> >& A, Matrix<Complex<R> >& t )
{
#ifndef RELEASE
    PushCallStack("HermitianTridiagU");
#endif
    const int tHeight = std::max(A.Height()-1,0);
#ifndef RELEASE
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( t.Viewing() && (t.Height() != tHeight || t.Width() != 1) )
        throw std::logic_error("t is of the wrong size");
#endif
    typedef Complex<R> C;
    
    if( !t.Viewing() )
        t.ResizeTo( tHeight, 1 );

    // Matrix views 
    Matrix<C>
        ATL, ATR,  A00, a01,     A02,  a01T,
        ABL, ABR,  a10, alpha11, a12,  alpha01B,
                   A20, a21,     A22;

    // Temporary matrices
    Matrix<C> w01;

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
        const C tau = Reflector( alpha01B, a01T );
        const R epsilon1 = alpha01B.GetRealPart(0,0);
        t.Set(t.Height()-1-A22.Height(),0,tau);
        alpha01B.Set(0,0,C(1));

        Hemv( UPPER, tau, A00, a01, C(0), w01 );
        const C alpha = -tau*Dot( w01, a01 )/C(2);
        Axpy( alpha, a01, w01 );
        Her2( UPPER, C(-1), a01, w01, A00 );
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

} // namespace hermitian_tridiag

template<typename R> 
inline void HermitianTridiag( UpperOrLower uplo, Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("HermitianTridiag");
#endif
    if( IsComplex<R>::val )
        throw std::logic_error("Called real routine with complex datatype");
    if( uplo == LOWER )
        hermitian_tridiag::HermitianTridiagL( A );
    else
        hermitian_tridiag::HermitianTridiagU( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void HermitianTridiag
( UpperOrLower uplo, Matrix<Complex<R> >& A, Matrix<Complex<R> >& t )
{
#ifndef RELEASE
    PushCallStack("HermitianTridiag");
#endif
    if( uplo == LOWER )
        hermitian_tridiag::HermitianTridiagL( A, t );
    else
        hermitian_tridiag::HermitianTridiagU( A, t );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
