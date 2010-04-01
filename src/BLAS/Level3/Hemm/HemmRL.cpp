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
#include "ElementalBLASInternal.h"
using namespace std;
using namespace Elemental;

template<typename T>
void
Elemental::BLAS::Internal::HemmRL
( const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::HemmRL");
#endif
    BLAS::Internal::HemmRLC( alpha, A, B, beta, C );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::BLAS::Internal::HemmRLC
( const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::HemmRLC");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "{A,B,C} must be distributed over the same grid.";
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        ATL(grid), ATR(grid),  A00(grid), A01(grid), A02(grid),  AColPan(grid),
        ABL(grid), ABR(grid),  A10(grid), A11(grid), A12(grid),  ARowPan(grid),
                               A20(grid), A21(grid), A22(grid);

    DistMatrix<T,MC,MR> BL(grid), BR(grid),
                        B0(grid), B1(grid), B2(grid);

    DistMatrix<T,MC,MR> CL(grid), CR(grid),
                        C0(grid), C1(grid), C2(grid),
                        CLeft(grid), CRight(grid);

    // Temporary distributions
    DistMatrix<T,MC,Star> B1_MC_Star(grid);
    DistMatrix<T,MR,Star> AColPan_MR_Star(grid);
    DistMatrix<T,Star,MR> ARowPan_Star_MR(grid);

    // Start the algorithm
    BLAS::Scal( beta, C );
    LockedPartitionDownDiagonal( A, ATL, ATR,
                                    ABL, ABR );
    LockedPartitionRight( B, BL, BR );
    PartitionRight( C, CL, CR );
    while( CR.Width() > 0 )
    {
        LockedRepartitionDownDiagonal( ATL, /**/ ATR,  A00, /**/ A01, A02,
                                      /*************/ /******************/
                                            /**/       A10, /**/ A11, A12,
                                       ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionRight( BL, /**/ BR,
                                B0, /**/ B1, B2 );

        RepartitionRight( CL, /**/ CR,
                          C0, /**/ C1, C2 );

        ARowPan.LockedView1x2( A10, A11 );

        AColPan.LockedView2x1( A11,
                               A21 );

        CLeft.View1x2( C0, C1 );

        CRight.View1x2( C1, C2 );

        B1_MC_Star.AlignWith( C );
        AColPan_MR_Star.AlignWith( CRight );
        ARowPan_Star_MR.AlignWith( CLeft );
        //--------------------------------------------------------------------//
        B1_MC_Star = B1;

        ARowPan_Star_MR = ARowPan;
        AColPan_MR_Star = AColPan;
        ARowPan_Star_MR.MakeTrapezoidal( Right, Lower );
        AColPan_MR_Star.MakeTrapezoidal( Left, Lower, -1 );

        BLAS::Gemm( Normal, Normal,
                    alpha, B1_MC_Star.LockedLocalMatrix(),
                           ARowPan_Star_MR.LockedLocalMatrix(),
                    (T)1,  CLeft.LocalMatrix()                 );
        BLAS::Gemm( Normal, ConjugateTranspose,
                    alpha, B1_MC_Star.LockedLocalMatrix(),
                           AColPan_MR_Star.LockedLocalMatrix(),
                    (T)1,  CRight.LocalMatrix()                );
        //--------------------------------------------------------------------//
        B1_MC_Star.FreeConstraints();
        AColPan_MR_Star.FreeConstraints();
        ARowPan_Star_MR.FreeConstraints();

        SlideLockedPartitionDownDiagonal( ATL, /**/ ATR,  A00, A01, /**/ A02,
                                               /**/       A10, A11, /**/ A12,
                                         /*************/ /******************/
                                          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionRight( BL,     /**/ BR,
                                   B0, B1, /**/ B2 );

        SlidePartitionRight( CL,     /**/ CR,
                             C0, C1, /**/ C2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Internal::HemmRL
( const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::HemmRL
( const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Internal::HemmRL
( const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                        const DistMatrix<scomplex,MC,MR>& B,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::HemmRL
( const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                        const DistMatrix<dcomplex,MC,MR>& B,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif
