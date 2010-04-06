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
#include "Elemental/BLASInternal.hpp"
using namespace std;
using namespace Elemental;

template<typename T>
void
Elemental::BLAS::Internal::HerkUC
( const T alpha, const DistMatrix<T,MC,MR>& A,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::HerkUC");
    if( A.GetGrid() != C.GetGrid() )
        throw "A and C must be distributed over the same grid.";
    if( A.Width() != C.Height() || A.Width() != C.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal HerkUT:" << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        const string s = msg.str();
        throw s.c_str();
    }
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> AT(grid),  A0(grid),
                        AB(grid),  A1(grid),
                                   A2(grid);

    // Temporary distributions
    DistMatrix<T,MR,  Star> A1_MR_Star(grid);
    DistMatrix<T,Star,VR  > A1Conj_Star_VR(grid);
    DistMatrix<T,Star,MC  > A1Conj_Star_MC(grid);

    // Start the algorithm
    BLAS::Scal( beta, C );
    LockedPartitionUp( A, AT, 
                          AB );
    while( AT.Height() > 0 )
    {
        LockedRepartitionUp( AT,  A0,
                                  A1,
                            /**/ /**/
                             AB,  A2 );

        A1_MR_Star.AlignWith( C );
        A1Conj_Star_MC.AlignWith( C );
        //--------------------------------------------------------------------//
        A1_MR_Star.TransposeFrom( A1 );
        A1Conj_Star_VR.ConjugateTransposeFrom( A1_MR_Star );
        A1Conj_Star_MC = A1Conj_Star_VR;

        BLAS::Internal::TriangularRankK
        ( Upper, alpha, A1Conj_Star_MC, A1_MR_Star, (T)1, C );
        //--------------------------------------------------------------------//
        A1_MR_Star.FreeConstraints();
        A1Conj_Star_MC.FreeConstraints();

        SlideLockedPartitionUp( AT,  A0,
                               /**/ /**/
                                     A1,
                                AB,  A2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Internal::HerkUC
( const float alpha, const DistMatrix<float,MC,MR>& A,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::HerkUC
( const double alpha, const DistMatrix<double,MC,MR>& A,
  const double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Internal::HerkUC
( const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::HerkUC
( const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

