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
Elemental::BLAS::Internal::Syr2kLN
( const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Syr2kLN");
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
        throw msg.str();
    }
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views 
    DistMatrix<T,MC,MR> AL(grid), AR(grid),
                        A0(grid), A1(grid), A2(grid);

    DistMatrix<T,MC,MR> BL(grid), BR(grid),
                        B0(grid), B1(grid), B2(grid);

    // Temporary distributions
    DistMatrix<T,MC,Star> A1_MC_Star(grid);
    DistMatrix<T,MR,Star> A1_MR_Star(grid);
    DistMatrix<T,MC,Star> B1_MC_Star(grid);
    DistMatrix<T,MR,Star> B1_MR_Star(grid);

    // Start the algorithm
    BLAS::Scal( beta, C );
    LockedPartitionRight( A, AL, AR );
    LockedPartitionRight( B, BL, BR );
    while( AR.Width() > 0 )
    {
        LockedRepartitionRight( AL, /**/ AR,
                                A0, /**/ A1, A2 );

        LockedRepartitionRight( BL, /**/ BR,
                                B0, /**/ B1, B2 );

        A1_MC_Star.AlignWith( C );
        A1_MR_Star.AlignWith( C );
        B1_MC_Star.AlignWith( C );
        B1_MR_Star.AlignWith( C );
        //--------------------------------------------------------------------//
        A1_MC_Star = A1;
        A1_MR_Star = A1_MC_Star;
        B1_MC_Star = B1;
        B1_MR_Star = B1_MC_Star;

        BLAS::Internal::Syr2kLNUpdate
        ( alpha, A1_MC_Star, A1_MR_Star, B1_MC_Star, B1_MR_Star, (T)1, C );
        //--------------------------------------------------------------------//
        A1_MC_Star.FreeConstraints();
        A1_MR_Star.FreeConstraints();
        B1_MC_Star.FreeConstraints();
        B1_MR_Star.FreeConstraints();

        SlideLockedPartitionRight( AL,     /**/ AR,
                                   A0, A1, /**/ A2 );

        SlideLockedPartitionRight( BL,     /**/ BR,
                                   B0, B1, /**/ B2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::BLAS::Internal::Syr2kLNUpdate
( const T alpha, const DistMatrix<T,MC,Star>& A_MC_Star,
                 const DistMatrix<T,MR,Star>& A_MR_Star,
                 const DistMatrix<T,MC,Star>& B_MC_Star,
                 const DistMatrix<T,MR,Star>& B_MR_Star,
  const T beta,        DistMatrix<T,MC,MR  >& C         )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Syr2kLNUpdate");
    if( A_MC_Star.GetGrid() != A_MR_Star.GetGrid() ||
        A_MR_Star.GetGrid() != B_MC_Star.GetGrid() ||
        B_MC_Star.GetGrid() != B_MR_Star.GetGrid() ||
        B_MR_Star.GetGrid() != C.GetGrid()           )
    {
        throw "{A,B,C} must be distributed over the same grid.";
    }
    if( A_MC_Star.Height() != C.Height() ||
        A_MR_Star.Height() != C.Width()  ||
        B_MC_Star.Height() != C.Height() ||
        B_MR_Star.Height() != C.Width()  ||
        A_MC_Star.Width()  != B_MR_Star.Width() ||
        B_MC_Star.Width()  != B_MR_Star.Width()    )
    {
        ostringstream msg;
        msg << "Nonconformal Syr2kLNUpdate: " << endl
            << "  A[MC,* ] ~ " << A_MC_Star.Height() << " x "
                               << A_MC_Star.Width()  << endl
            << "  A[MR,* ] ~ " << A_MR_Star.Height() << " x "
                               << A_MR_Star.Width()  << endl
            << "  B[MC,* ] ~ " << B_MC_Star.Height() << " x " 
                               << B_MC_Star.Width()  << endl
            << "  B[MR,* ] ~ " << B_MR_Star.Height() << " x "
                               << B_MR_Star.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw msg.str();
    }
    if( A_MC_Star.ColAlignment() != C.ColAlignment() ||
        A_MR_Star.ColAlignment() != C.RowAlignment() ||
        B_MC_Star.ColAlignment() != C.ColAlignment() ||
        B_MR_Star.ColAlignment() != C.RowAlignment()    )
    {
        ostringstream msg;
        msg << "Misaligned Syr2kLNUpdate: " << endl
            << "  A[MC,* ] ~ " << A_MC_Star.ColAlignment() << endl
            << "  A[MR,* ] ~ " << A_MR_Star.ColAlignment() << endl
            << "  B[MC,* ] ~ " << B_MC_Star.ColAlignment() << endl
            << "  B[MR,* ] ~ " << B_MR_Star.ColAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw msg.str();
    }
#endif
    const Grid& grid = C.GetGrid();

    // Matrix views 
    DistMatrix<T,MC,Star> AT_MC_Star(grid),  A0_MC_Star(grid),
                          AB_MC_Star(grid),  A1_MC_Star(grid),
                                             A2_MC_Star(grid);

    DistMatrix<T,MR,Star> AT_MR_Star(grid),  A0_MR_Star(grid),
                          AB_MR_Star(grid),  A1_MR_Star(grid),
                                             A2_MR_Star(grid);

    DistMatrix<T,MC,Star> BT_MC_Star(grid),  B0_MC_Star(grid),
                          BB_MC_Star(grid),  B1_MC_Star(grid),
                                             B2_MC_Star(grid);

    DistMatrix<T,MR,Star> BT_MR_Star(grid),  B0_MR_Star(grid),
                          BB_MR_Star(grid),  B1_MR_Star(grid),
                                             B2_MR_Star(grid);

    DistMatrix<T,MC,MR> 
        CTL(grid), CTR(grid),  C00(grid), C01(grid), C02(grid),
        CBL(grid), CBR(grid),  C10(grid), C11(grid), C12(grid),
                               C20(grid), C21(grid), C22(grid);

    DistMatrix<T,MC,MR> D11(grid);

    // We want our local gemms to be of width blocksize, and so we will 
    // temporarily change to c times the current blocksize
    PushBlocksizeStack( grid.Width()*Blocksize() );

    // Start the algorithm
    BLAS::Scal( beta, C );
    LockedPartitionDown( A_MC_Star, AT_MC_Star,
                                    AB_MC_Star );
    LockedPartitionDown( A_MR_Star, AT_MR_Star,
                                    AB_MR_Star );
    LockedPartitionDown( B_MC_Star, BT_MC_Star,
                                    BB_MC_Star );
    LockedPartitionDown( B_MR_Star, BT_MR_Star,
                                    BB_MR_Star );
    PartitionDownDiagonal( C, CTL, CTR,
                              CBL, CBR );
    while( AB_MC_Star.Height() > 0 )
    {
        LockedRepartitionDown( AT_MC_Star,  A0_MC_Star,
                              /**********/ /**********/
                                            A1_MC_Star,
                               AB_MC_Star,  A2_MC_Star );

        LockedRepartitionDown( AT_MR_Star,  A0_MR_Star,
                              /**********/ /**********/
                                            A1_MR_Star,
                               AB_MR_Star,  A2_MR_Star );

        LockedRepartitionDown( BT_MC_Star,  B0_MC_Star,
                              /**********/ /**********/
                                            B1_MC_Star,
                               BB_MC_Star,  B2_MC_Star );

        LockedRepartitionDown( BT_MR_Star,  B0_MR_Star,
                              /**********/ /**********/
                                            B1_MR_Star,
                               BB_MR_Star,  B2_MR_Star );

        RepartitionDownDiagonal( CTL, /**/ CTR,  C00, /**/ C01, C02,
                                /*************/ /******************/
                                      /**/       C10, /**/ C11, C12,
                                 CBL, /**/ CBR,  C20, /**/ C21, C22 );

        D11.AlignWith( C11 );
        D11.ResizeTo( C11.Height(), C11.Width() );
        //--------------------------------------------------------------------//
        BLAS::Gemm( Normal, Transpose,
                    alpha, A1_MC_Star.LockedLocalMatrix(),
                           B1_MR_Star.LockedLocalMatrix(),
                    (T)0,  D11.LocalMatrix()              );
        BLAS::Gemm( Normal, Transpose,
                    alpha, A2_MC_Star.LockedLocalMatrix(),
                           B1_MR_Star.LockedLocalMatrix(),
                    (T)1,  C21.LocalMatrix()              );

        BLAS::Gemm( Normal, Transpose,
                    alpha, B2_MC_Star.LockedLocalMatrix(),
                           A1_MR_Star.LockedLocalMatrix(),
                    (T)1,  C21.LocalMatrix()              );
        BLAS::Gemm( Normal, Transpose,
                    alpha, B1_MC_Star.LockedLocalMatrix(),
                           A1_MR_Star.LockedLocalMatrix(),
                    (T)1,  D11.LocalMatrix()              );

        D11.MakeTrapezoidal( Left, Lower );
        BLAS::Axpy( (T)1, D11, C11 );
        //--------------------------------------------------------------------//
        D11.FreeConstraints();

        SlideLockedPartitionDown( AT_MC_Star,  A0_MC_Star,
                                               A1_MC_Star,
                                 /**********/ /**********/
                                  AB_MC_Star,  A2_MC_Star );

        SlideLockedPartitionDown( AT_MR_Star,  A0_MR_Star,
                                               A1_MR_Star,
                                 /**********/ /**********/
                                  AB_MR_Star,  A2_MR_Star );

        SlideLockedPartitionDown( BT_MC_Star,  B0_MC_Star,
                                               B1_MC_Star,
                                 /**********/ /**********/
                                  BB_MC_Star,  B2_MC_Star );

        SlideLockedPartitionDown( BT_MR_Star,  B0_MR_Star,
                                               B1_MR_Star,
                                 /**********/ /**********/
                                  BB_MR_Star,  B2_MR_Star );

        SlidePartitionDownDiagonal( CTL, /**/ CTR,  C00, C01, /**/ C02,
                                         /**/       C10, C11, /**/ C12,
                                   /*************/ /******************/
                                    CBL, /**/ CBR,  C20, C21, /**/ C22 );
    }

    PopBlocksizeStack();

#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Internal::Syr2kLN
( const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::Syr2kLN
( const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Internal::Syr2kLN
( const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                        const DistMatrix<scomplex,MC,MR>& B,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::Syr2kLN
( const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                        const DistMatrix<dcomplex,MC,MR>& B,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

