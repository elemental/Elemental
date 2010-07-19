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
elemental::blas::internal::HerkLC
( T alpha, const DistMatrix<T,MC,MR>& A,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::HerkLC");
    if( A.GetGrid() != C.GetGrid() )
        throw logic_error( "A and C must be distributed over the same grid." );
    if( A.Width() != C.Height() || A.Width() != C.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal HerkLC:" << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> AT(g),  A0(g),
                        AB(g),  A1(g),
                                A2(g);

    // Temporary distributions
    DistMatrix<T,MR,  Star> A1Trans_MR_Star(g);
    DistMatrix<T,Star,VR  > A1_Star_VR(g);
    DistMatrix<T,Star,MC  > A1_Star_MC(g);

    // Start the algorithm
    blas::Scal( beta, C );
    LockedPartitionDown
    ( A, AT, 
         AB, 0 );
    while( AB.Height() > 0 )
    {
        LockedRepartitionDown
        ( AT,  A0,
         /**/ /**/
               A1,
          AB,  A2 );

        A1Trans_MR_Star.AlignWith( C );
        A1_Star_MC.AlignWith( C );
        //--------------------------------------------------------------------//
        A1Trans_MR_Star.TransposeFrom( A1 );
        A1_Star_VR.TransposeFrom( A1Trans_MR_Star );
        A1_Star_MC = A1_Star_VR;

        blas::internal::LocalTriangularRankK
        ( Lower, ConjugateTranspose, Transpose,
          alpha, A1_Star_MC, A1Trans_MR_Star, (T)1, C );
        //--------------------------------------------------------------------//
        A1Trans_MR_Star.FreeAlignments();
        A1_Star_MC.FreeAlignments();

        SlideLockedPartitionDown
        ( AT,  A0,
               A1,
         /**/ /**/
          AB,  A2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::HerkLC
( float alpha, const DistMatrix<float,MC,MR>& A,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::HerkLC
( double alpha, const DistMatrix<double,MC,MR>& A,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::HerkLC
( scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::HerkLC
( dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

