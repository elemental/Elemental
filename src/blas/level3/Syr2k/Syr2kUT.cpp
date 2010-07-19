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
elemental::blas::internal::Syr2kUT
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::Syr2kUT");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
    if( A.Width() != C.Height() || 
        A.Width() != C.Width()  ||
        B.Width() != C.Height() ||
        B.Width() != C.Width()  ||
        A.Height() != B.Height()  )
    {
        ostringstream msg;
        msg << "Nonconformal Syr2kUT:" << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B ~ " << B.Height() << " x " << B.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> AT(g),  A0(g),
                        AB(g),  A1(g),
                                A2(g);

    DistMatrix<T,MC,MR> BT(g),  B0(g),
                        BB(g),  B1(g),
                                B2(g);

    // Temporary distributions
    DistMatrix<T,MR,  Star> A1Trans_MR_Star(g);
    DistMatrix<T,MR,  Star> B1Trans_MR_Star(g);
    DistMatrix<T,Star,VR  > A1_Star_VR(g);
    DistMatrix<T,Star,VR  > B1_Star_VR(g);
    DistMatrix<T,Star,MC  > A1_Star_MC(g);
    DistMatrix<T,Star,MC  > B1_Star_MC(g);

    // Start the algorithm
    blas::Scal( beta, C );
    LockedPartitionDown
    ( A, AT, 
         AB, 0 );
    LockedPartitionDown
    ( B, BT,
         BB, 0 );
    while( AB.Height() > 0 )
    {
        LockedRepartitionDown
        ( AT,  A0,
         /**/ /**/
               A1,
          AB,  A2 );

        LockedRepartitionDown
        ( BT,  B0,
         /**/ /**/
               B1,
          BB,  B2 );

        A1Trans_MR_Star.AlignWith( C );
        B1Trans_MR_Star.AlignWith( C );
        A1_Star_MC.AlignWith( C );
        B1_Star_MC.AlignWith( C );
        //--------------------------------------------------------------------//
        A1Trans_MR_Star.TransposeFrom( A1 );
        A1_Star_VR.TransposeFrom( A1Trans_MR_Star );
        A1_Star_MC = A1_Star_VR;

        B1Trans_MR_Star.TransposeFrom( B1 );
        B1_Star_VR.TransposeFrom( B1Trans_MR_Star );
        B1_Star_MC = B1_Star_VR;

        blas::internal::LocalTriangularRank2K
        ( Upper, Transpose, Transpose, Transpose, Transpose, alpha,
          A1_Star_MC, B1_Star_MC, A1Trans_MR_Star, B1Trans_MR_Star, (T)1, C );
        //--------------------------------------------------------------------//
        A1Trans_MR_Star.FreeAlignments();
        B1Trans_MR_Star.FreeAlignments();
        A1_Star_MC.FreeAlignments();
        B1_Star_MC.FreeAlignments();

        SlideLockedPartitionDown
        ( AT,  A0,
               A1,
         /**/ /**/
          AB,  A2 );

        SlideLockedPartitionDown
        ( BT,  B0,
               B1,
         /**/ /**/
          BB,  B2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::Syr2kUT
( float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::Syr2kUT
( double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::Syr2kUT
( scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::Syr2kUT
( dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

