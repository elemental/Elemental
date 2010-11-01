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

#ifdef TIMING
namespace {
bool alreadyTimed = false;

Timer timer;
Timer getA21_VC_Star_Timer;
Timer getA11_Star_Star_Timer;
Timer panelLUTimer;
Timer composePivotsTimer;
Timer applyPivotsTimer;
Timer localTrsmTimer;
Timer getA21Trans_Star_MC_Timer;
Timer getA12_Star_VR_Timer;
Timer localGemmTimer;
Timer storeResultsTimer;

inline void
ResetTimers()
{
    ::timer.Reset();
    ::getA21_VC_Star_Timer.Reset();
    ::getA11_Star_Star_Timer.Reset();
    ::panelLUTimer.Reset();
    ::composePivotsTimer.Reset();
    ::applyPivotsTimer.Reset();
    ::localTrsmTimer.Reset();
    ::getA21Trans_Star_MC_Timer.Reset();
    ::getA12_Star_VR_Timer.Reset();
    ::localGemmTimer.Reset();
    ::storeResultsTimer.Reset();
}
} // anonymous namespace

void
elemental::lapack::lu::PrintTimings()
{
#ifndef RELEASE
    PushCallStack("lapack::lu::PrintTimings");
    if( !::alreadyTimed )
        throw std::logic_error("You have not yet run elemental::lapack::LU");
#endif
    cout << "\n"
         << "elemental::lapack::LU breakdown:\n" 
         << "--------------------------------------------------------\n"
         << "A21[VC,* ] <- A21[MC,MR]:     "
         << ::getA21_VC_Star_Timer.Time() << " seconds.\n"
         << "A11[* ,* ] <- A11[MC,MR]:     "
         << ::getA11_Star_Star_Timer.Time() << " seconds.\n"
         << "Panel LU:                     "
         << ::panelLUTimer.Time() << " seconds.\n"
         << "Pivot composition:            "
         << ::composePivotsTimer.Time() << " seconds.\n"
         << "Pivot application:            "
         << ::applyPivotsTimer.Time() << " seconds.\n"
         << "Local Trsm:                   "
         << ::localTrsmTimer.Time() << " seconds.\n"
         << "(A21^T)[* ,MC] <- A21[VC,* ]: "
         << ::getA21Trans_Star_MC_Timer.Time() << " seconds.\n"
         << "A12[* ,VR] <- A12[MC,MR]:     "
         << ::getA12_Star_VR_Timer.Time() << " seconds.\n"
         << "Local Gemm:                   "
         << ::localGemmTimer.Time() << " seconds.\n"
         << "Storing results:              "
         << ::storeResultsTimer.Time() << " seconds.\n"
         << "Total: " << ::timer.Time() << " seconds.\n" << std::endl;
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // TIMING

template<typename T>
void
elemental::lapack::LU
( DistMatrix<T,MC,MR>& A, DistMatrix<int,VC,Star>& p )
{
#ifndef RELEASE
    PushCallStack("lapack::LU");
    if( A.GetGrid() != p.GetGrid() )
        throw logic_error( "A and p must be distributed over the same grid." );
    if( A.Height() != p.Height() ) 
        throw logic_error( "A and p must be the same height." );
#endif
#ifdef TIMING
    ::ResetTimers();
    ::timer.Start();
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  AB(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),  
                         A20(g), A21(g), A22(g);

    DistMatrix<int,VC,Star>
        pT(g),  p0(g), 
        pB(g),  p1(g),
                p2(g);

    // Temporary distributions
    DistMatrix<T,  Star,Star> A11_Star_Star(g);
    DistMatrix<T,  Star,VR  > A12_Star_VR(g);
    DistMatrix<T,  Star,MR  > A12_Star_MR(g);
    DistMatrix<T,  VC,  Star> A21_VC_Star(g);
    DistMatrix<T,  Star,MC  > A21Trans_Star_MC(g);
    DistMatrix<int,Star,Star> p1_Star_Star(g);

    // Pivot composition
    vector<int> image;
    vector<int> preimage;

    // Start the algorithm
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( p, pT,
         pB, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        RepartitionDown
        ( pT,  p0,
         /**/ /**/
               p1,
          pB,  p2 );

        AB.View1x2( ABL, ABR );

        int pivotOffset = A01.Height();
        A12_Star_VR.AlignWith( A22 );
        A12_Star_MR.AlignWith( A22 );
        A21_VC_Star.AlignWith( A22 );
        A21Trans_Star_MC.AlignWith( A22 );
        A11_Star_Star.ResizeTo( A11.Height(), A11.Width() );
        p1_Star_Star.ResizeTo( p1.Height(), 1 );
        //--------------------------------------------------------------------//
#ifdef TIMING
        ::getA21_VC_Star_Timer.Start();
#endif
        A21_VC_Star = A21;
#ifdef TIMING
        ::getA21_VC_Star_Timer.Stop();
        ::getA11_Star_Star_Timer.Start();
#endif
        A11_Star_Star = A11;
#ifdef TIMING
        ::getA11_Star_Star_Timer.Stop();
        ::panelLUTimer.Start();
#endif
        lapack::internal::PanelLU
        ( A11_Star_Star, 
          A21_VC_Star, p1_Star_Star, pivotOffset );
#ifdef TIMING
        ::panelLUTimer.Stop();
        ::composePivotsTimer.Start();
#endif
        lapack::internal::ComposePivots
        ( p1_Star_Star, image, preimage, pivotOffset );
#ifdef TIMING
        ::composePivotsTimer.Stop();
        ::applyPivotsTimer.Start();
#endif
        lapack::internal::ApplyRowPivots( AB, image, preimage, pivotOffset );
#ifdef TIMING
        ::applyPivotsTimer.Stop();
        ::getA12_Star_VR_Timer.Start();
#endif
        A12_Star_VR = A12;
#ifdef TIMING
        ::getA12_Star_VR_Timer.Stop();
        ::localTrsmTimer.Start();
#endif
        blas::internal::LocalTrsm
        ( Left, Lower, Normal, Unit, (T)1, A11_Star_Star, A12_Star_VR );
#ifdef TIMING
        ::localTrsmTimer.Stop();
        ::getA21Trans_Star_MC_Timer.Start();
#endif
        A21Trans_Star_MC.TransposeFrom( A21_VC_Star );
#ifdef TIMING
        ::getA21Trans_Star_MC_Timer.Stop();
        ::getA12_Star_VR_Timer.Start();
#endif
        A12_Star_MR = A12_Star_VR;
#ifdef TIMING
        ::getA12_Star_VR_Timer.Stop();
        ::localGemmTimer.Start();
#endif
        blas::internal::LocalGemm
        ( Transpose, Normal, (T)-1, A21Trans_Star_MC, A12_Star_MR, (T)1, A22 );
#ifdef TIMING
        ::localGemmTimer.Stop();
        ::storeResultsTimer.Start();
#endif
        A11 = A11_Star_Star;
        A12 = A12_Star_MR;
        A21.TransposeFrom( A21Trans_Star_MC );
        p1 = p1_Star_Star;
#ifdef TIMING
        ::storeResultsTimer.Stop();
#endif
        //--------------------------------------------------------------------//
        A12_Star_VR.FreeAlignments();
        A12_Star_MR.FreeAlignments();
        A21_VC_Star.FreeAlignments();
        A21Trans_Star_MC.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlidePartitionDown
        ( pT,  p0,
               p1,
         /**/ /**/
          pB,  p2 );
    }
#ifdef TIMING
    ::timer.Stop();
    ::alreadyTimed = true;
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template void
elemental::lapack::LU
( DistMatrix<float,MC,MR>& A, DistMatrix<int,VC,Star>& p );

template void
elemental::lapack::LU
( DistMatrix<double,MC,MR>& A, DistMatrix<int,VC,Star>& p );

#ifndef WITHOUT_COMPLEX
template void
elemental::lapack::LU
( DistMatrix<scomplex,MC,MR>& A, DistMatrix<int,VC,Star>& p );

template void
elemental::lapack::LU
( DistMatrix<dcomplex,MC,MR>& A, DistMatrix<int,VC,Star>& p );
#endif

