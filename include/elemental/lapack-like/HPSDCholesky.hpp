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

#ifndef WITHOUT_PMRRR

namespace elem {
namespace hpsd_cholesky {

template<typename F>
void MakeExplicitlyHermitian( UpperOrLower uplo, DistMatrix<F>& A )
{
    const Grid& g = A.Grid();
    DistMatrix<F> ATL(g), ATR(g),  A00(g), A01(g), A02(g),
                  ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                                   A20(g), A21(g), A22(g);
    DistMatrix<F> A11Adj(g);
    DistMatrix<F,MR,MC> A11_MR_MC(g);
    DistMatrix<F,MR,MC> A21_MR_MC(g);
    DistMatrix<F,MR,MC> A12_MR_MC(g);

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        A11Adj.AlignWith( A11 );
        A11_MR_MC.AlignWith( A11 );
        A12_MR_MC.AlignWith( A21 );
        A21_MR_MC.AlignWith( A12 );
        //--------------------------------------------------------------------//
        A11_MR_MC = A11;
        A11Adj.ResizeTo( A11.Height(), A11.Width() );
        Adjoint( A11_MR_MC.LocalMatrix(), A11Adj.LocalMatrix() );

        if( uplo == LOWER )
        {
            MakeTrapezoidal( LEFT, UPPER, 1, A11Adj );
            Axpy( F(1), A11Adj, A11 );

            A21_MR_MC = A21;
            Adjoint( A21_MR_MC.LocalMatrix(), A12.LocalMatrix() ); 
        }
        else
        {
            MakeTrapezoidal( LEFT, LOWER, -1, A11Adj );
            Axpy( F(1), A11Adj, A11 );

            A12_MR_MC = A12;
            Adjoint( A12_MR_MC.LocalMatrix(), A21.LocalMatrix() );
        }
        //--------------------------------------------------------------------//
        A21_MR_MC.FreeAlignments();
        A12_MR_MC.FreeAlignments();
        A11_MR_MC.FreeAlignments();
        A11Adj.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
}

} // namespace hpsd_cholesky

//
// Compute the Cholesky factor of a potentially singular Hermitian semi-definite
// matrix.
//

template<typename R>
inline void
HPSDCholesky( UpperOrLower uplo, DistMatrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("HPSDCholesky");
#endif
    HPSDSquareRoot( uplo, A );
    hpsd_cholesky::MakeExplicitlyHermitian( uplo, A );

    if( uplo == LOWER )
    {
        LQ( A );
        MakeTrapezoidal( LEFT, LOWER, 0, A );
    }
    else
    {
        QR( A );
        MakeTrapezoidal( RIGHT, UPPER, 0, A );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
HPSDCholesky( UpperOrLower uplo, DistMatrix<Complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("HPSDCholesky");
#endif
    HPSDSquareRoot( uplo, A );
    hpsd_cholesky::MakeExplicitlyHermitian( uplo, A );

    const Grid& g = A.Grid();
    if( uplo == LOWER )
    {
        DistMatrix<Complex<R>,MD,STAR> t(g);
        LQ( A, t );
        MakeTrapezoidal( LEFT, LOWER, 0, A );
    }
    else
    {
        DistMatrix<Complex<R>,MD,STAR> t(g);
        QR( A, t );
        MakeTrapezoidal( RIGHT, UPPER, 0, A );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // WITHOUT_PMRRR
