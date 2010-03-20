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
#include "ElementalBLAS_Internal.h"
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
#endif
    const Grid& grid = A.GetGrid();
#ifndef RELEASE
    if( A.GetGrid() != C.GetGrid() )
    {
        if( grid.VCRank() == 0 )
            cerr << "A and C must be distributed over the same grid." << endl;
        DumpCallStack();
        throw exception();
    }
    if( A.Width() != C.Height() || A.Width() != C.Width() )
    {
        if( grid.VCRank() == 0 )
        {
            cerr << "Nonconformal HerkUT:" <<
            endl << "  A ~ " << A.Height() << " x " << A.Width() <<
            endl << "  C ~ " << C.Height() << " x " << C.Width() << 
            endl;
        }
        DumpCallStack();
        throw exception();
    }
#endif
    // Matrix views
    DistMatrix<T,MC,MR> AT(grid),  A0(grid),
                        AB(grid),  A1(grid),
                                   A2(grid);

    // Temporary distributions
    DistMatrix<T,Star,MC> A1_Star_MC(grid);
    DistMatrix<T,Star,MR> A1_Star_MR(grid);

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

        A1_Star_MC.AlignWith( C );
        A1_Star_MR.AlignWith( C );
        //--------------------------------------------------------------------//
        A1_Star_MR = A1;
        A1_Star_MC = A1_Star_MR;

        BLAS::Internal::HerkUC_Update
        ( alpha, A1_Star_MC, A1_Star_MR, (T)1, C );
        //--------------------------------------------------------------------//
        A1_Star_MC.FreeConstraints();
        A1_Star_MR.FreeConstraints();

        SlideLockedPartitionUp( AT,  A0,
                               /**/ /**/
                                     A1,
                                AB,  A2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::BLAS::Internal::HerkUC_Update
( const T alpha, const DistMatrix<T,Star,MC>& A_Star_MC,
                 const DistMatrix<T,Star,MR>& A_Star_MR,
  const T beta,        DistMatrix<T,MC,  MR>& C         )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::HerkUC_Update");
#endif
    const Grid& grid = C.GetGrid();
#ifndef RELEASE
    if( A_Star_MC.GetGrid() != A_Star_MR.GetGrid() || 
        A_Star_MR.GetGrid() != C.GetGrid()            )
    {
        if( grid.VCRank() == 0 )
            cerr << "A and C must be distributed over the same grid." << endl;
        DumpCallStack();
        throw exception();
    }
    if( A_Star_MC.Width() != C.Height() ||
        A_Star_MR.Width() != C.Width()  ||
        A_Star_MC.Height() != A_Star_MR.Height() ||
        A_Star_MC.Width()  != A_Star_MR.Width()    )
    {
        if( grid.VCRank() == 0 )
        {
            cerr << "Nonconformal HerkUC_Update: " <<
            endl << "  A[* ,MC] ~ " << A_Star_MC.Height() << " x "
                                    << A_Star_MC.Width()  <<
            endl << "  A[* ,MR] ~ " << A_Star_MR.Height() << " x "
                                    << A_Star_MR.Width()  <<
            endl << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        }
        DumpCallStack();
        throw exception();
    }
    if( A_Star_MC.RowAlignment() != C.ColAlignment() ||
        A_Star_MR.RowAlignment() != C.RowAlignment()   )
    {
        if( grid.VCRank() == 0 )
        {
            cerr << "Misaligned HerkUC_Update: " <<
            endl << "  A[* ,MC] ~ " << A_Star_MC.RowAlignment() <<
            endl << "  A[* ,MR] ~ " << A_Star_MR.RowAlignment() <<
            endl << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                       C.RowAlignment() << endl;
        }
        DumpCallStack();
        throw exception();
    }
#endif
    // Matrix views
    DistMatrix<T,Star,MC> 
        AL_Star_MC(grid), AR_Star_MC(grid),
        A0_Star_MC(grid), A1_Star_MC(grid), A2_Star_MC(grid);

    DistMatrix<T,Star,MR> 
        AL_Star_MR(grid), AR_Star_MR(grid),
        A0_Star_MR(grid), A1_Star_MR(grid), A2_Star_MR(grid);

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
    LockedPartitionLeft( A_Star_MC, AL_Star_MC, AR_Star_MC );
    LockedPartitionLeft( A_Star_MR, AL_Star_MR, AR_Star_MR );
    PartitionUpDiagonal( C, CTL, CTR,
                            CBL, CBR );
    while( AL_Star_MC.Width() > 0 )
    {
        LockedRepartitionLeft( AL_Star_MC,             /**/ AR_Star_MC,
                               A0_Star_MC, A1_Star_MC, /**/ A2_Star_MC );

        LockedRepartitionLeft( AL_Star_MR,             /**/ AR_Star_MR,
                               A0_Star_MR, A1_Star_MR, /**/ A2_Star_MR );

        RepartitionUpDiagonal( CTL, /**/ CTR,  C00, C01, /**/ C02,
                                    /**/       C10, C11, /**/ C12,
                              /*************/ /******************/
                               CBL, /**/ CBR,  C20, C21, /**/ C22 );

        D11.AlignWith( C11 );
        D11.ResizeTo( C11.Height(), C11.Width() );
        //--------------------------------------------------------------------//
        BLAS::Gemm( ConjugateTranspose, Normal, 
                    alpha, A0_Star_MC.LockedLocalMatrix(),
                           A1_Star_MR.LockedLocalMatrix(),
                    (T)1,  C01.LocalMatrix()              );
        BLAS::Gemm( ConjugateTranspose, Normal,
                    alpha, A1_Star_MC.LockedLocalMatrix(),
                           A1_Star_MR.LockedLocalMatrix(),
                    (T)0,  D11.LocalMatrix()              );
        D11.MakeTrapezoidal( Left, Upper );
        BLAS::Axpy( (T)1, D11, C11 );
        //--------------------------------------------------------------------//
        D11.FreeConstraints();
        
        SlideLockedPartitionLeft( AL_Star_MC, /**/ AR_Star_MC,
                                  A0_Star_MC, /**/ A1_Star_MC, A2_Star_MC );

        SlideLockedPartitionLeft( AL_Star_MR, /**/ AR_Star_MR,
                                  A0_Star_MR, /**/ A1_Star_MR, A2_Star_MR );

        SlidePartitionUpDiagonal( CTL, /**/ CTR,  C00, /**/ C01, C02,
                                 /*************/ /******************/
                                       /**/       C10, /**/ C11, C12,
                                  CBL, /**/ CBR,  C20, /**/ C21, C22 );
    }

    PopBlocksizeStack();

#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Internal::HerkUC
( const float alpha, const DistMatrix<float,MC,MR>& A,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::HerkUC_Update
( const float alpha, const DistMatrix<float,Star,MC>& A_Star_MC,
                     const DistMatrix<float,Star,MR>& A_Star_MR,
  const float beta,        DistMatrix<float,MC,  MR>& C         );

template void Elemental::BLAS::Internal::HerkUC
( const double alpha, const DistMatrix<double,MC,MR>& A,
  const double beta,        DistMatrix<double,MC,MR>& C );

template void Elemental::BLAS::Internal::HerkUC_Update
( const double alpha, const DistMatrix<double,Star,MC>& A_Star_MC,
                      const DistMatrix<double,Star,MR>& A_Star_MR,
  const double beta,        DistMatrix<double,MC,  MR>& C         );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Internal::HerkUC
( const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::HerkUC_Update
( const scomplex alpha, 
  const DistMatrix<scomplex,Star,MC>& A_Star_MC,
  const DistMatrix<scomplex,Star,MR>& A_Star_MR,
  const scomplex beta,        
        DistMatrix<scomplex,MC,  MR>& C         );

template void Elemental::BLAS::Internal::HerkUC
( const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::HerkUC_Update
( const dcomplex alpha, 
  const DistMatrix<dcomplex,Star,MC>& A_Star_MC,
  const DistMatrix<dcomplex,Star,MR>& A_Star_MR,
  const dcomplex beta,        
        DistMatrix<dcomplex,MC,  MR>& C         );
#endif

