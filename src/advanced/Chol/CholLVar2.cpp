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

/*
   Parallelization of Variant 2 Lower Cholesky factorization. 

   Original serial update:
   ------------------------
   A11 := A11 - A10 A10^H
   A11 := Chol(A11)
   A21 := A21 - A20 A10^H
   A21 := A21 tril(A11)^-H
   ------------------------

   Our parallel update:
   -----------------------------------------------------
   A10^H[MR,* ] <- A10[MC,MR]
   X11[MC,* ] := A10[MC,MR] A10^H[MR,* ]
   A11[MC,MR] := A11[MC,MR] - (SumRow(X11[MC,* ]))[* ,MR]

   A11[* ,* ] <- A11[MC,MR]   
   A11[* ,* ] := Chol(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]

   X21[MC,* ] := A20[MC,MR] A10^H[MR,* ]
   A21[MC,MR] := A21[MC,MR] - (SumRow(X21[MC,* ]))[* ,MR]

   A21[VC,* ] <- A21[MC,MR]
   A21[VC,* ] := A21[VC,* ] tril(A11[* ,* ])^-H
   A21[MC,MR] <- A21[VC,* ]
   -----------------------------------------------------
*/
template<typename F> // representation of real or complex number
void
elemental::advanced::internal::CholLVar2
( DistMatrix<F,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::CholLVar2");
    if( A.Height() != A.Width() )
        throw logic_error
        ( "Can only compute Cholesky factor of square matrices." );
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F,MC,MR> 
        ATL(g), ATR(g),   A00(g), A01(g), A02(g),
        ABL(g), ABR(g),   A10(g), A11(g), A12(g),
                          A20(g), A21(g), A22(g);

    // Temporary distributions
    DistMatrix<F,MR,  Star> A10Herm_MR_Star(g);
    DistMatrix<F,Star,Star> A11_Star_Star(g);
    DistMatrix<F,VC,  Star> A21_VC_Star(g);
    DistMatrix<F,MC,  Star> X11_MC_Star(g);
    DistMatrix<F,MC,  Star> X21_MC_Star(g);

    // Start the algorithm
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

        A10Herm_MR_Star.AlignWith( A10 );
        X11_MC_Star.AlignWith( A10 );
        X21_MC_Star.AlignWith( A20 );
        X11_MC_Star.ResizeTo( A11.Height(), A11.Width() );
        X21_MC_Star.ResizeTo( A21.Height(), A21.Width() );
        //--------------------------------------------------------------------//
        A10Herm_MR_Star.ConjugateTransposeFrom( A10 );
        basic::internal::LocalGemm
        ( Normal, Normal, 
          (F)1, A10, A10Herm_MR_Star, (F)0, X11_MC_Star );
        A11.SumScatterUpdate( (F)-1, X11_MC_Star );

        A11_Star_Star = A11;
        advanced::internal::LocalChol( Lower, A11_Star_Star );
        A11 = A11_Star_Star;

        basic::internal::LocalGemm
        ( Normal, Normal,
          (F)1, A20, A10Herm_MR_Star, (F)0, X21_MC_Star );
        A21.SumScatterUpdate( (F)-1, X21_MC_Star );

        A21_VC_Star = A21;
        basic::internal::LocalTrsm
        ( Right, Lower, ConjugateTranspose, NonUnit, 
          (F)1, A11_Star_Star, A21_VC_Star );
        A21 = A21_VC_Star;
        //--------------------------------------------------------------------//
        A10Herm_MR_Star.FreeAlignments();
        X11_MC_Star.FreeAlignments();
        X21_MC_Star.FreeAlignments();

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

/*
   Naive parallelization of Variant 2 Lower Cholesky factorization. 

   Original serial update:
   ------------------------
   A11 := A11 - A10 A10^H
   A11 := Chol(A11)
   A21 := A21 - A20 A10^H
   A21 := A21 tril(A11)^-H
   ------------------------

   Our parallel update:
   -----------------------------------------------------
   A10[* ,MR] <- A10[MC,MR]
   X11[MC,* ] := A10[MC,MR] (A10[* ,MR])^H
   A11[MC,MR] := A11[MC,MR] - (SumRow(X11[MC,* ]))[* ,MR]

   A11[* ,* ] <- A11[MC,MR]   
   A11[* ,* ] := Chol(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]

   X21[MC,* ] := A20[MC,MR] (A10[* ,MR])^H
   A21[MC,MR] := A21[MC,MR] - (SumRow(X21[MC,* ]))[* ,MR]

   A21[VC,* ] <- A21[MC,MR]
   A21[VC,* ] := A21[VC,* ] tril(A11[* ,* ])^-H
   A21[MC,MR] <- A21[VC,* ]
   -----------------------------------------------------
*/
template<typename F> // representation of real or complex number
void
elemental::advanced::internal::CholLVar2Naive
( DistMatrix<F,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::CholLVar2Naive");
    if( A.Height() != A.Width() )
        throw logic_error
        ( "Can only compute Cholesky factor of square matrices." );
    if( A.Grid().VCRank() == 0 )
    {
        cout << "CholLVar2Naive exists solely for academic purposes. Please "
                "use CholLVar2 in real applications." << endl;
    }
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F,MC,MR> 
        ATL(g), ATR(g),   A00(g), A01(g), A02(g),
        ABL(g), ABR(g),   A10(g), A11(g), A12(g),
                          A20(g), A21(g), A22(g);

    // Temporary distributions
    DistMatrix<F,Star,MR  > A10_Star_MR(g);
    DistMatrix<F,Star,Star> A11_Star_Star(g);
    DistMatrix<F,VC,  Star> A21_VC_Star(g);
    DistMatrix<F,MC,  Star> X11_MC_Star(g);
    DistMatrix<F,MC,  Star> X21_MC_Star(g);

    // Start the algorithm
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

        A10_Star_MR.AlignWith( A10 );
        X11_MC_Star.AlignWith( A10 );
        X21_MC_Star.AlignWith( A20 );
        X11_MC_Star.ResizeTo( A11.Height(), A11.Width() );
        X21_MC_Star.ResizeTo( A21.Height(), A21.Width() );
        //--------------------------------------------------------------------//
        A10_Star_MR = A10;
        basic::internal::LocalGemm
        ( Normal, ConjugateTranspose, 
          (F)1, A10, A10_Star_MR, (F)0, X11_MC_Star );
        A11.SumScatterUpdate( (F)-1, X11_MC_Star );

        A11_Star_Star = A11;
        advanced::internal::LocalChol( Lower, A11_Star_Star );
        A11 = A11_Star_Star;

        basic::internal::LocalGemm
        ( Normal, ConjugateTranspose, 
          (F)1, A20, A10_Star_MR, (F)0, X21_MC_Star );
        A21.SumScatterUpdate( (F)-1, X21_MC_Star );

        A21_VC_Star = A21;
        basic::internal::LocalTrsm
        ( Right, Lower, ConjugateTranspose, NonUnit, 
          (F)1, A11_Star_Star, A21_VC_Star );
        A21 = A21_VC_Star;
        //--------------------------------------------------------------------//
        A10_Star_MR.FreeAlignments();
        X11_MC_Star.FreeAlignments();
        X21_MC_Star.FreeAlignments();

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

template void elemental::advanced::internal::CholLVar2
( DistMatrix<float,MC,MR>& A );

template void elemental::advanced::internal::CholLVar2Naive
( DistMatrix<float,MC,MR>& A );

template void elemental::advanced::internal::CholLVar2
( DistMatrix<double,MC,MR>& A );

template void elemental::advanced::internal::CholLVar2Naive
( DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::advanced::internal::CholLVar2
( DistMatrix<scomplex,MC,MR>& A );

template void elemental::advanced::internal::CholLVar2Naive
( DistMatrix<scomplex,MC,MR>& A );

template void elemental::advanced::internal::CholLVar2
( DistMatrix<dcomplex,MC,MR>& A );

template void elemental::advanced::internal::CholLVar2Naive
( DistMatrix<dcomplex,MC,MR>& A );
#endif

