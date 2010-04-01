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
Elemental::BLAS::Internal::SyrkUT
( const T alpha, const DistMatrix<T,MC,MR>& A,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::SyrkUT");
    if( A.GetGrid() != C.GetGrid() )
        throw "A and C must be distributed over the same grid.";
    if( A.Width() != C.Height() || A.Width() != C.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal SyrkUT:" << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        throw msg.str();
    }
#endif
    const Grid& grid = A.GetGrid();

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

        BLAS::Internal::SyrkUTUpdate
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
Elemental::BLAS::Internal::SyrkUTUpdate
( const T alpha, const DistMatrix<T,Star,MC>& A_Star_MC,
                 const DistMatrix<T,Star,MR>& A_Star_MR,
  const T beta,        DistMatrix<T,MC,  MR>& C         )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::SyrkUTUpdate");
    if( A_Star_MC.GetGrid() != A_Star_MR.GetGrid() || 
        A_Star_MR.GetGrid() != C.GetGrid()           )
    {
        throw "A and C must be distributed over the same grid.";
    }
    if( A_Star_MC.Width() != C.Height() ||
        A_Star_MR.Width() != C.Width()  ||
        A_Star_MC.Height() != A_Star_MR.Height() ||
        A_Star_MC.Width()  != A_Star_MR.Width()    )
    {
        ostringstream msg;
        msg << "Nonconformal SyrkUTUpdate: " << endl
            << "  A[* ,MC] ~ " << A_Star_MC.Height() << " x "
                               << A_Star_MC.Width()  << endl
            << "  A[* ,MR] ~ " << A_Star_MR.Height() << " x "
                               << A_Star_MR.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw msg.str();
    }
    if( A_Star_MC.RowAlignment() != C.ColAlignment() ||
        A_Star_MR.RowAlignment() != C.RowAlignment()   )
    {
        ostringstream msg;
        msg << "Misaligned SyrkUTUpdate: " << endl
            << "  A[* ,MC] ~ " << A_Star_MC.RowAlignment() << endl
            << "  A[* ,MR] ~ " << A_Star_MR.RowAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw msg.str();
    }
#endif
    const Grid& grid = C.GetGrid();

    if( C.Height() < 2*grid.Width()*Blocksize() )
    {
        BLAS::Internal::SyrkUTUpdateKernel
        ( alpha, A_Star_MC, A_Star_MR, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on CTR
        // and recurse on CTL and CBR.

        DistMatrix<T,Star,MC> AL_Star_MC(grid), AR_Star_MC(grid);
        DistMatrix<T,Star,MR> AL_Star_MR(grid), AR_Star_MR(grid);
        DistMatrix<T,MC,MR> CTL(grid), CTR(grid),
                            CBL(grid), CBR(grid);

        const unsigned half = C.Height() / 2;

        LockedPartitionRight( A_Star_MC, AL_Star_MC, AR_Star_MC, half );
        LockedPartitionRight( A_Star_MR, AL_Star_MR, AR_Star_MR, half );
        PartitionDownDiagonal( C, CTL, CTR,
                                  CBL, CBR, half );

        BLAS::Gemm
        ( Transpose, Normal,
          alpha, AL_Star_MC.LockedLocalMatrix(),
                 AR_Star_MR.LockedLocalMatrix(),
          beta,  CTR.LocalMatrix()              );

        // Recurse
        BLAS::Internal::SyrkUTUpdate
        ( alpha, AL_Star_MC, AL_Star_MR, beta, CTL );

        BLAS::Internal::SyrkUTUpdate
        ( alpha, AR_Star_MC, AR_Star_MR, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::BLAS::Internal::SyrkUTUpdateKernel
( const T alpha, const DistMatrix<T,Star,MC>& A_Star_MC,
                 const DistMatrix<T,Star,MR>& A_Star_MR,
  const T beta,        DistMatrix<T,MC,  MR>& C         )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::SyrkUTUpdateKernel");
    if( A_Star_MC.GetGrid() != A_Star_MR.GetGrid() || 
        A_Star_MR.GetGrid() != C.GetGrid()           )
    {
        throw "A and C must be distributed over the same grid.";
    }
    if( A_Star_MC.Width() != C.Height() ||
        A_Star_MR.Width() != C.Width()  ||
        A_Star_MC.Height() != A_Star_MR.Height() ||
        A_Star_MC.Width()  != A_Star_MR.Width()    )
    {
        ostringstream msg;
        msg << "Nonconformal SyrkUTUpdateKernel: " << endl
            << "  A[* ,MC] ~ " << A_Star_MC.Height() << " x "
                               << A_Star_MC.Width()  << endl
            << "  A[* ,MR] ~ " << A_Star_MR.Height() << " x "
                               << A_Star_MR.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw msg.str();
    }
    if( A_Star_MC.RowAlignment() != C.ColAlignment() ||
        A_Star_MR.RowAlignment() != C.RowAlignment()   )
    {
        ostringstream msg;
        msg << "Misaligned SyrkUTUpdateKernel: " << endl
            << "  A[* ,MC] ~ " << A_Star_MC.RowAlignment() << endl
            << "  A[* ,MR] ~ " << A_Star_MR.RowAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw msg.str();
    }
#endif
    const Grid& grid = C.GetGrid();

    DistMatrix<T,Star,MC> AL_Star_MC(grid), AR_Star_MC(grid);
    DistMatrix<T,Star,MR> AL_Star_MR(grid), AR_Star_MR(grid);

    DistMatrix<T,MC,MR>
        CTL(grid), CTR(grid),
        CBL(grid), CBR(grid);

    DistMatrix<T,MC,MR> DTL(grid), DBR(grid);

    const unsigned half = C.Height()/2;

    BLAS::Scal( beta, C );

    LockedPartitionRight( A_Star_MC, AL_Star_MC, AR_Star_MC, half );
    LockedPartitionRight( A_Star_MR, AL_Star_MR, AR_Star_MR, half );
    PartitionDownDiagonal( C, CTL, CTR,
                              CBL, CBR, half );

    DTL.AlignWith( CTL );
    DBR.AlignWith( CBR );
    DTL.ResizeTo( CTL.Height(), CTL.Width() );
    DBR.ResizeTo( CBR.Height(), CBR.Width() );
    //------------------------------------------------------------------------//
    BLAS::Gemm( Transpose, Normal,
                alpha, AL_Star_MC.LockedLocalMatrix(),
                       AR_Star_MR.LockedLocalMatrix(),
                (T)1,  CTR.LocalMatrix()              );

    BLAS::Gemm( Transpose, Normal,
                alpha, AL_Star_MC.LockedLocalMatrix(),
                       AL_Star_MR.LockedLocalMatrix(),
                (T)0,  DTL.LocalMatrix()              );
    DTL.MakeTrapezoidal( Left, Upper );
    BLAS::Axpy( (T)1, DTL, CTL );

    BLAS::Gemm( Transpose, Normal,
                alpha, AR_Star_MC.LockedLocalMatrix(),
                       AR_Star_MR.LockedLocalMatrix(),
                (T)0,  DBR.LocalMatrix()              );
    DBR.MakeTrapezoidal( Left, Upper );
    BLAS::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//

#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Internal::SyrkUT
( const float alpha, const DistMatrix<float,MC,MR>& A,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::SyrkUT
( const double alpha, const DistMatrix<double,MC,MR>& A,
  const double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Internal::SyrkUT
( const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::SyrkUT
( const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

