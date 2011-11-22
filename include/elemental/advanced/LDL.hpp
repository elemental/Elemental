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

#include "./LDL/LocalLDL.hpp"

template<typename F>
inline void elemental::advanced::LDLT
( DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,STAR>& d )
{
#ifndef RELEASE
    PushCallStack("advanced::LDLT");
#endif
    advanced::internal::LDLVar3( TRANSPOSE, A, d );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void elemental::advanced::LDLH
( DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,STAR>& d )
{
#ifndef RELEASE
    PushCallStack("advanced::LDLH");
#endif
    advanced::internal::LDLVar3( ADJOINT, A, d );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> // F represents a real or complex field
inline void
elemental::advanced::internal::LDLVar3
( Orientation orientation, DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,STAR>& d )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::LDLVar3");
    if( orientation == NORMAL )
        throw std::logic_error("Can only perform LDL^T and LDL^H");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Grid() != d.Grid() )
        throw std::logic_error("A and d must use the same grid");
    if( d.Viewing() && (d.Height() != A.Height() || d.Width() != 1) )
        throw std::logic_error
        ("d must be a column vector of the same height as A");
    if( d.Viewing() && d.ColAlignment() != A.ColAlignment() )
        throw std::logic_error("d must be aligned with A");
#endif
    const Grid& g = A.Grid();
    if( !d.Viewing() )
    {
        d.AlignWith( A );
        d.ResizeTo( A.Height(), 1 );
    }

    // Matrix views
    DistMatrix<F,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<F,MC,STAR>
        dT(g),  d0(g),
        dB(g),  d1(g),
                d2(g);

    // Temporary matrices
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> d1_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> A21_VR_STAR(g);
    DistMatrix<F,STAR,MC  > S21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A21AdjOrTrans_STAR_MR(g);

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( d, dT,
         dB, 0 );
    while( ABR.Height() > 0 )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        RepartitionDown
        ( dT,  d0,
         /**/ /**/
               d1,
          dB,  d2 );

        A21_VC_STAR.AlignWith( A22 );
        A21_VR_STAR.AlignWith( A22 );
        S21Trans_STAR_MC.AlignWith( A22 );
        A21AdjOrTrans_STAR_MR.AlignWith( A22 );
        //--------------------------------------------------------------------//
        A11_STAR_STAR = A11; 
        advanced::internal::LocalLDL
        ( orientation, A11_STAR_STAR, d1_STAR_STAR );
        A11 = A11_STAR_STAR;
        d1 = d1_STAR_STAR;

        A21_VC_STAR = A21;
        basic::internal::LocalTrsm
        ( RIGHT, LOWER, orientation, UNIT, 
          (F)1, A11_STAR_STAR, A21_VC_STAR );

        S21Trans_STAR_MC.TransposeFrom( A21_VC_STAR );
        basic::DiagonalSolve( RIGHT, NORMAL, d1_STAR_STAR, A21_VC_STAR );
        A21_VR_STAR = A21_VC_STAR;
        if( orientation == ADJOINT )
            A21AdjOrTrans_STAR_MR.AdjointFrom( A21_VR_STAR );
        else
            A21AdjOrTrans_STAR_MR.TransposeFrom( A21_VR_STAR );

        basic::internal::LocalTrrk
        ( LOWER, TRANSPOSE,
          (F)-1, S21Trans_STAR_MC, A21AdjOrTrans_STAR_MR, (F)1, A22 );

        basic::DiagonalSolve( LEFT, NORMAL, d1_STAR_STAR, S21Trans_STAR_MC );
        A21.TransposeFrom( S21Trans_STAR_MC );
        //--------------------------------------------------------------------//
        A21_VC_STAR.FreeAlignments();
        A21_VR_STAR.FreeAlignments();
        S21Trans_STAR_MC.FreeAlignments();
        A21AdjOrTrans_STAR_MR.FreeAlignments();

        SlidePartitionDown
        ( dT,  d0,
               d1,
         /**/ /**/
          dB,  d2 );

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
