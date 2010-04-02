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
Elemental::BLAS::Internal::SymmLL
( const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::SymmLL");
#endif
    BLAS::Internal::SymmLLC( alpha, A, B, beta, C );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::BLAS::Internal::SymmLLC
( const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::SymmLLC");
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
    BLAS::Scal( beta, C );
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

        ARowPan.LockedView1x2( A10, A11 );

        AColPan.LockedView2x1( A11,
                               A21 );

        CAbove.View2x1( C0,
                        C1 );

        CBelow.View2x1( C1,
                        C2 );

        AColPan_MC_Star.AlignWith( CBelow );
        ARowPan_Star_MC.AlignWith( CAbove );
        B1_Star_MR.AlignWith( C );
        //--------------------------------------------------------------------//
        AColPan_MC_Star = AColPan;
        ARowPan_Star_MC = ARowPan;
        AColPan_MC_Star.MakeTrapezoidal( Left, Lower );
        ARowPan_Star_MC.MakeTrapezoidal( Right, Lower, -1 );

        B1_Star_MR = B1;

        BLAS::Gemm( Normal, Normal,
                    alpha, AColPan_MC_Star.LockedLocalMatrix(),
                           B1_Star_MR.LockedLocalMatrix(),
                    (T)1,  CBelow.LocalMatrix()                );
        BLAS::Gemm( Transpose, Normal,
                    alpha, ARowPan_Star_MC.LockedLocalMatrix(),
                           B1_Star_MR.LockedLocalMatrix(),
                    (T)1,  CAbove.LocalMatrix()                );
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

template void Elemental::BLAS::Internal::SymmLL
( const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::SymmLL
( const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Internal::SymmLL
( const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                        const DistMatrix<scomplex,MC,MR>& B,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::SymmLL
( const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                        const DistMatrix<dcomplex,MC,MR>& B,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

