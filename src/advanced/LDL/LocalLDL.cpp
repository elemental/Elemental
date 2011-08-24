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
#include "elemental/basic_internal.hpp"
#include "elemental/advanced_internal.hpp"
using namespace std;
using namespace elemental;

template<typename F>
void
elemental::advanced::LDLH( Matrix<F>& A, Matrix<F>& d )
{
#ifndef RELEASE
    PushCallStack("advanced::LDLH");
#endif
    advanced::internal::LDLVar3( ADJOINT, A, d );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
void
elemental::advanced::LDLT( Matrix<F>& A, Matrix<F>& d )
{
#ifndef RELEASE
    PushCallStack("advanced::LDLT");
#endif
    advanced::internal::LDLVar3( TRANSPOSE, A, d );
#ifndef RELEASE
    PopCallStack();
#endif
}

// Unblocked serial LDL _without_ partial pivoting
template<typename F> // representation of a field
void
elemental::advanced::internal::LDLVar3
( Orientation orientation, Matrix<F>& A, Matrix<F>& d )
{
#ifndef RELEASE
    PushCallStack("advanced::LDLVar3");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( d.Viewing() && (d.Height() != A.Height() || d.Width() != 1) )
        throw std::logic_error
        ("d must be a column vector the same height as A");
    if( orientation == NORMAL )
        throw std::logic_error("Can only perform LDL^T or LDL^H");
#endif
    if( !d.Viewing() )
        d.ResizeTo( A.Height(), 1 );

    // Matrix views 
    Matrix<F>
        ATL, ATR,  A00, a01,     A02,  alpha21T,
        ABL, ABR,  a10, alpha11, a12,  a21B,
                   A20, a21,     A22;
    Matrix<F>
        dT,  d0,
        dB,  delta1,
             d2;

    // Temporary matrices
    Matrix<F> s21;

    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( d, dT,
         dB, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        RepartitionDown
        ( dT,  d0,
         /**/ /******/
               delta1,
          dB,  d2 );

        //--------------------------------------------------------------------//
        F delta = alpha11.Get(0,0);
        delta1.Set(0,0,delta);

        s21 = a21;
        basic::Scal( static_cast<F>(1)/delta, a21 );
        if( orientation == ADJOINT )
            basic::Ger( (F)-1, a21, s21, A22 );
        else
            basic::Geru( (F)-1, a21, s21, A22 );
        //--------------------------------------------------------------------//

        SlidePartitionDown
        ( dT,  d0,
               delta1,
         /**/ /******/
          dB,  d2 );

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

template void elemental::advanced::LDLH( Matrix<float>& A, Matrix<float>& d );
template void elemental::advanced::LDLT( Matrix<float>& A, Matrix<float>& d );
template void elemental::advanced::internal::LDLVar3
( Orientation orientation, Matrix<float>& A, Matrix<float>& d );

template void elemental::advanced::LDLH( Matrix<double>& A, Matrix<double>& d );
template void elemental::advanced::LDLT( Matrix<double>& A, Matrix<double>& d );
template void elemental::advanced::internal::LDLVar3
( Orientation orientation, Matrix<double>& A, Matrix<double>& d );

#ifndef WITHOUT_COMPLEX
template void elemental::advanced::LDLH
( Matrix<scomplex>& A, Matrix<scomplex>& d );
template void elemental::advanced::LDLT
( Matrix<scomplex>& A, Matrix<scomplex>& d );
template void elemental::advanced::internal::LDLVar3
( Orientation orientation, Matrix<scomplex>& A, Matrix<scomplex>& d );

template void elemental::advanced::LDLH
( Matrix<dcomplex>& A, Matrix<dcomplex>& d );
template void elemental::advanced::LDLT
( Matrix<dcomplex>& A, Matrix<dcomplex>& d );
template void elemental::advanced::internal::LDLVar3
( Orientation orientation, Matrix<dcomplex>& A, Matrix<dcomplex>& d );
#endif
