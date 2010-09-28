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
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::blas::internal::SyrkLN
( T alpha, const DistMatrix<T,MC,MR>& A,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::SyrkLN");
    if( A.GetGrid() != C.GetGrid() )
        throw logic_error( "A and C must be distributed over the same grid." );
    if( A.Height() != C.Height() || A.Height() != C.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal SyrkLN:" << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width() <<  endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> AL(g), AR(g),
                        A0(g), A1(g), A2(g);

    // Temporary distributions
    DistMatrix<T,MC,  Star> A1_MC_Star(g);
    DistMatrix<T,VR,  Star> A1_VR_Star(g);
    DistMatrix<T,Star,MR  > A1Trans_Star_MR(g);

    // Start the algorithm
    C.ScaleTrapezoidal( beta, Left, Lower );
    LockedPartitionRight( A, AL, AR, 0 );
    while( AR.Width() > 0 )
    {
        LockedRepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );

        A1_MC_Star.AlignWith( C );
        A1_VR_Star.AlignWith( C );
        A1Trans_Star_MR.AlignWith( C );
        //--------------------------------------------------------------------//
        A1_VR_Star = A1_MC_Star = A1;
        A1Trans_Star_MR.TransposeFrom( A1_VR_Star );

        blas::internal::LocalTriangularRankK
        ( Lower, alpha, A1_MC_Star, A1Trans_Star_MR, (T)1, C );
        //--------------------------------------------------------------------//
        A1_MC_Star.FreeAlignments();
        A1_VR_Star.FreeAlignments();
        A1Trans_Star_MR.FreeAlignments();

        SlideLockedPartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::SyrkLN
( float alpha, const DistMatrix<float,MC,MR>& A,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::SyrkLN
( double alpha, const DistMatrix<double,MC,MR>& A,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::SyrkLN
( scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::SyrkLN
( dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

