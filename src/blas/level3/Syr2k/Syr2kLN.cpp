/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::blas::internal::Syr2kLN
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::Syr2kLN");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "{A,B,C} must be distributed over the same grid.";
    if( A.Height() != C.Height() || A.Height() != C.Width() ||
        B.Height() != C.Height() || B.Height() != C.Width() ||
        A.Width() != B.Width()                                 )
    {
        ostringstream msg;
        msg << "Nonconformal Syr2kLN:" << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B ~ " << B.Height() << " x " << B.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        const string& s = msg.str();
        throw s.c_str();
    }
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views 
    DistMatrix<T,MC,MR> AL(grid), AR(grid),
                        A0(grid), A1(grid), A2(grid);

    DistMatrix<T,MC,MR> BL(grid), BR(grid),
                        B0(grid), B1(grid), B2(grid);

    // Temporary distributions
    DistMatrix<T,MC,  Star> A1_MC_Star(grid);
    DistMatrix<T,MC,  Star> B1_MC_Star(grid);
    DistMatrix<T,VR,  Star> A1_VR_Star(grid);
    DistMatrix<T,VR,  Star> B1_VR_Star(grid);
    DistMatrix<T,Star,MR  > A1Trans_Star_MR(grid);
    DistMatrix<T,Star,MR  > B1Trans_Star_MR(grid);

    // Start the algorithm
    blas::Scal( beta, C );
    LockedPartitionRight( A, AL, AR );
    LockedPartitionRight( B, BL, BR );
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
        A1Trans_Star_MR.AlignWith( C );
        B1Trans_Star_MR.AlignWith( C );
        //--------------------------------------------------------------------//
        A1_VR_Star = A1_MC_Star = A1;
        A1Trans_Star_MR.TransposeFrom( A1_VR_Star );

        B1_VR_Star = B1_MC_Star = B1;
        B1Trans_Star_MR.TransposeFrom( B1_VR_Star );

        blas::internal::TriangularRank2K
        ( Lower, alpha, 
          A1_MC_Star, B1_MC_Star, A1Trans_Star_MR, B1Trans_Star_MR, (T)1, C );
        //--------------------------------------------------------------------//
        A1_MC_Star.FreeAlignments();
        B1_MC_Star.FreeAlignments();
        A1_VR_Star.FreeAlignments();
        B1_VR_Star.FreeAlignments();
        A1Trans_Star_MR.FreeAlignments();
        B1Trans_Star_MR.FreeAlignments();

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

template void elemental::blas::internal::Syr2kLN
( float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::Syr2kLN
( double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::Syr2kLN
( scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::Syr2kLN
( dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

