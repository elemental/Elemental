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
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::blas::internal::Her2kLN
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::Her2kLN");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
    if( A.Height() != C.Height() || A.Height() != C.Width() ||
        B.Height() != C.Height() || B.Height() != C.Width() ||
        A.Width() != B.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal Her2kLN:" << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B ~ " << B.Height() << " x " << B.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = A.GetGrid();

    // Matrix views 
    DistMatrix<T,MC,MR> AL(g), AR(g),
                        A0(g), A1(g), A2(g);

    DistMatrix<T,MC,MR> BL(g), BR(g),
                        B0(g), B1(g), B2(g);

    // Temporary distributions
    DistMatrix<T,MC,  Star> A1_MC_Star(g);
    DistMatrix<T,MC,  Star> B1_MC_Star(g);
    DistMatrix<T,VR,  Star> A1_VR_Star(g);
    DistMatrix<T,VR,  Star> B1_VR_Star(g);
    DistMatrix<T,Star,MR  > A1Herm_Star_MR(g);
    DistMatrix<T,Star,MR  > B1Herm_Star_MR(g);

    // Start the algorithm
    blas::Scal( beta, C );
    LockedPartitionRight( A, AL, AR, 0 );
    LockedPartitionRight( B, BL, BR, 0 );
    while( AR.Width() > 0 )
    {
        LockedRepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );

        LockedRepartitionRight
        ( BL, /**/ BR,
          B0, /**/ B1, B2 );

        A1_MC_Star.AlignWith( C );
        B1_MC_Star.AlignWith( C );
        A1_VR_Star.AlignWith( C );
        B1_VR_Star.AlignWith( C );
        A1Herm_Star_MR.AlignWith( C );
        B1Herm_Star_MR.AlignWith( C );
        //--------------------------------------------------------------------//
        A1_VR_Star = A1_MC_Star = A1;
        A1Herm_Star_MR.ConjugateTransposeFrom( A1_VR_Star );

        B1_VR_Star = B1_MC_Star = B1;
        B1Herm_Star_MR.ConjugateTransposeFrom( B1_VR_Star );

        blas::internal::LocalTriangularRank2K
        ( Lower, alpha,
          A1_MC_Star, B1_MC_Star, A1Herm_Star_MR, B1Herm_Star_MR, (T)1, C );
        //--------------------------------------------------------------------//
        A1_MC_Star.FreeAlignments();
        B1_MC_Star.FreeAlignments();
        A1_VR_Star.FreeAlignments();
        B1_VR_Star.FreeAlignments();
        A1Herm_Star_MR.FreeAlignments();
        B1Herm_Star_MR.FreeAlignments();

        SlideLockedPartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );

        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::Her2kLN
( float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::Her2kLN
( double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::Her2kLN
( scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::Her2kLN
( dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

