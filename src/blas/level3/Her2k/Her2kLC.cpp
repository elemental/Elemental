/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::blas::internal::Her2kLC
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::Her2kLC");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "{A,B,C} must be distributed over the same grid.";
    if( A.Width() != C.Height() || 
        A.Width() != C.Width()  ||
        B.Width() != C.Height() ||
        B.Width() != C.Width()  ||
        A.Height() != B.Height()  )
    {
        ostringstream msg;
        msg << "Nonconformal Her2kLC:" << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B ~ " << B.Height() << " x " << B.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        const string& s = msg.str();
        throw s.c_str();
    }
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> AT(grid),  A0(grid),
                        AB(grid),  A1(grid),
                                   A2(grid);

    DistMatrix<T,MC,MR> BT(grid),  B0(grid),
                        BB(grid),  B1(grid),
                                   B2(grid);

    // Temporary distributions
    DistMatrix<T,MR,  Star> A1_MR_Star(grid);
    DistMatrix<T,MR,  Star> B1_MR_Star(grid);
    DistMatrix<T,Star,VR  > A1Conj_Star_VR(grid);
    DistMatrix<T,Star,VR  > B1Conj_Star_VR(grid);
    DistMatrix<T,Star,MC  > A1Conj_Star_MC(grid);
    DistMatrix<T,Star,MC  > B1Conj_Star_MC(grid);

    // Start the algorithm
    blas::Scal( beta, C );
    LockedPartitionDown( A, AT,
                            AB );
    LockedPartitionDown( B, BT,
                            BB );
    while( AB.Height() > 0 )
    {
        LockedRepartitionDown( AT,  A0,
                              /**/ /**/
                                    A1,
                               AB,  A2 );

        LockedRepartitionDown( BT,  B0,
                              /**/ /**/
                                    B1,
                               BB,  B2 );

        A1_MR_Star.AlignWith( C );
        B1_MR_Star.AlignWith( C );
        A1Conj_Star_MC.AlignWith( C );
        B1Conj_Star_MC.AlignWith( C );
        //--------------------------------------------------------------------//
        A1_MR_Star.TransposeFrom( A1 );
        A1Conj_Star_VR.ConjugateTransposeFrom( A1_MR_Star );
        A1Conj_Star_MC = A1Conj_Star_VR;

        B1_MR_Star.TransposeFrom( B1 );
        B1Conj_Star_VR.ConjugateTransposeFrom( B1_MR_Star );
        B1Conj_Star_MC = B1Conj_Star_VR;

        blas::internal::TriangularRank2K
        ( Lower, alpha,
          A1Conj_Star_MC, B1Conj_Star_MC, A1_MR_Star, B1_MR_Star, (T)1, C );
        //--------------------------------------------------------------------//
        A1_MR_Star.FreeConstraints();
        B1_MR_Star.FreeConstraints();
        A1Conj_Star_MC.FreeConstraints();
        B1Conj_Star_MC.FreeConstraints();

        SlideLockedPartitionDown( AT,  A0,
                                       A1,
                                 /**/ /**/
                                  AB,  A2 );

        SlideLockedPartitionDown( BT,  B0,
                                       B1,
                                 /**/ /**/
                                  BB,  B2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::Her2kLC
( float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::Her2kLC
( double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::Her2kLC
( scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::Her2kLC
( dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

