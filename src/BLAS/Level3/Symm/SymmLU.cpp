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
elemental::blas::internal::SymmLU
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::SymmLU");
#endif
    blas::internal::SymmLUC( alpha, A, B, beta, C );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::SymmLUC
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::SymmLUC");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "{A,B,C} must be distributed over the same grid.";
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        ATL(grid), ATR(grid),  A00(grid), A01(grid), A02(grid),  AColPan(grid),
        ABL(grid), ABR(grid),  A10(grid), A11(grid), A12(grid),  ARowPan(grid),
                               A20(grid), A21(grid), A22(grid);

    DistMatrix<T,MC,MR> BT(grid),  B0(grid),
                        BB(grid),  B1(grid),
                                   B2(grid);

    DistMatrix<T,MC,MR> CT(grid),  C0(grid),  CAbove(grid),
                        CB(grid),  C1(grid),  CBelow(grid),
                                   C2(grid);

    // Temporary distributions
    DistMatrix<T,MC,Star> AColPan_MC_Star(grid);
    DistMatrix<T,Star,MC> ARowPan_Star_MC(grid);
    DistMatrix<T,Star,MR> B1_Star_MR(grid);

    // Start the algorithm
    blas::Scal( beta, C );
    LockedPartitionDownDiagonal( A, ATL, ATR,
                                    ABL, ABR );
    LockedPartitionDown( B, BT,
                            BB );
    PartitionDown( C, CT,
                      CB );
    while( CB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal( ATL, /**/ ATR,  A00, /**/ A01, A02,
                                      /*************/ /******************/
                                            /**/       A10, /**/ A11, A12,
                                       ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionDown( BT,  B0,
                              /**/ /**/
                                    B1,
                               BB,  B2 );

        RepartitionDown( CT,  C0,
                        /**/ /**/
                              C1,
                         CB,  C2 );

        ARowPan.LockedView1x2( A11, A12 );

        AColPan.LockedView2x1( A01,
                               A11 );

        CAbove.View2x1( C0,
                        C1 );

        CBelow.View2x1( C1,
                        C2 );

        AColPan_MC_Star.AlignWith( CAbove );
        ARowPan_Star_MC.AlignWith( CBelow );
        B1_Star_MR.AlignWith( C );
        //--------------------------------------------------------------------//
        AColPan_MC_Star = AColPan;
        ARowPan_Star_MC = ARowPan;
        AColPan_MC_Star.MakeTrapezoidal( Right, Upper );
        ARowPan_Star_MC.MakeTrapezoidal( Left, Upper, 1 );

        B1_Star_MR = B1;

        blas::Gemm( Normal, Normal,
                    alpha, AColPan_MC_Star.LockedLocalMatrix(),
                           B1_Star_MR.LockedLocalMatrix(),
                    (T)1,  CAbove.LocalMatrix()                );
        blas::Gemm( Transpose, Normal,
                    alpha, ARowPan_Star_MC.LockedLocalMatrix(),
                           B1_Star_MR.LockedLocalMatrix(),
                    (T)1,  CBelow.LocalMatrix()                );
        //--------------------------------------------------------------------//
        AColPan_MC_Star.FreeConstraints();
        ARowPan_Star_MC.FreeConstraints();
        B1_Star_MR.FreeConstraints();

        SlideLockedPartitionDownDiagonal( ATL, /**/ ATR,  A00, A01, /**/ A02,
                                               /**/       A10, A11, /**/ A12,
                                         /*************/ /******************/
                                          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDown( BT,  B0,
                                       B1,
                                 /**/ /**/
                                  BB,  B2 );

        SlidePartitionDown( CT,  C0,
                                 C1,
                           /**/ /**/
                            CB,  C2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::SymmLU
( float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::SymmLU
( double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::SymmLU
( scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::SymmLU
( dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

