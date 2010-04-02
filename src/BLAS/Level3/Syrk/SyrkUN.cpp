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
Elemental::BLAS::Internal::SyrkUN
( const T alpha, const DistMatrix<T,MC,MR>& A,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::SyrkUN");
    if( A.GetGrid() != C.GetGrid() )
        throw "A and C must be distributed over the same grid.";
    if( A.Height() != C.Height() || A.Height() != C.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal SyrkUN:" << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        throw msg.str();
    }
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> AL(grid), AR(grid),
                        A0(grid), A1(grid), A2(grid);

    // Temporary distributions
    DistMatrix<T,MC,Star> A1_MC_Star(grid);
    DistMatrix<T,MR,Star> A1_MR_Star(grid);

    // Start the algorithm
    BLAS::Scal( beta, C );
    LockedPartitionRight( A, AL, AR );
    while( AR.Width() > 0 )
    {
        LockedRepartitionRight( AL, /**/ AR,
                                A0, /**/ A1, A2 );

        A1_MC_Star.AlignWith( C );
        A1_MR_Star.AlignWith( C );
        //--------------------------------------------------------------------//
        A1_MC_Star = A1;
        A1_MR_Star = A1_MC_Star;

        BLAS::Internal::SyrkUNUpdate
        ( alpha, A1_MC_Star, A1_MR_Star, (T)1, C ); 
        //--------------------------------------------------------------------//
        A1_MC_Star.FreeConstraints();
        A1_MR_Star.FreeConstraints();

        SlideLockedPartitionRight( AL,     /**/ AR,
                                   A0, A1, /**/ A2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::BLAS::Internal::SyrkUNUpdate
( const T alpha, const DistMatrix<T,MC,Star>& A_MC_Star,
                 const DistMatrix<T,MR,Star>& A_MR_Star,
  const T beta,        DistMatrix<T,MC,MR  >& C         )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::SyrkUNUpdate");
    if( A_MC_Star.GetGrid() != A_MR_Star.GetGrid() || 
        A_MR_Star.GetGrid() != C.GetGrid()           )
    {
        throw "A and C must be distributed over the same grid.";
    }
    if( A_MC_Star.Height() != C.Height() ||
        A_MR_Star.Height() != C.Width()  ||
        A_MC_Star.Height() != A_MR_Star.Height() ||
        A_MC_Star.Width()  != A_MR_Star.Width()    )
    {
        ostringstream msg;
        msg << "Nonconformal SyrkUNUpdate: " << endl
            << "  A[MC,* ] ~ " << A_MC_Star.Height() << " x "
                               << A_MC_Star.Width()  << endl
            << "  A[MR,* ] ~ " << A_MR_Star.Height() << " x "
                               << A_MR_Star.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw msg.str();
    }
    if( A_MC_Star.ColAlignment() != C.ColAlignment() ||
        A_MR_Star.ColAlignment() != C.RowAlignment()   )
    {
        ostringstream msg;
        msg << "Misaligned SyrkUNUpdate: " << endl
            << "  A[MC,* ] ~ " << A_MC_Star.ColAlignment() << endl
            << "  A[MR,* ] ~ " << A_MR_Star.ColAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw msg.str();
    }
#endif
    const Grid& grid = C.GetGrid();

    if( C.Height() < 2*grid.Width()*Blocksize() )
    {
        BLAS::Internal::SyrkUNUpdateKernel
        ( alpha, A_MC_Star, A_MR_Star, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on CTR
        // and recurse on CTL and CBR.

        DistMatrix<T,MC,Star> AT_MC_Star(grid),
                              AB_MC_Star(grid);

        DistMatrix<T,MR,Star> AT_MR_Star(grid),
                              AB_MR_Star(grid);

        DistMatrix<T,MC,MR> CTL(grid), CTR(grid),
                            CBL(grid), CBR(grid);

        const unsigned half = C.Height() / 2;

        LockedPartitionDown( A_MC_Star, AT_MC_Star,
                                        AB_MC_Star, half );

        LockedPartitionDown( A_MR_Star, AT_MR_Star,
                                        AB_MR_Star, half );

        PartitionDownDiagonal( C, CTL, CTR,
                                  CBL, CBR, half );

        BLAS::Gemm
        ( Normal, Transpose,
          alpha, AT_MC_Star.LockedLocalMatrix(),
                 AB_MR_Star.LockedLocalMatrix(),
          beta,  CTR.LocalMatrix()              );

        // Recurse
        BLAS::Internal::SyrkUNUpdate
        ( alpha, AT_MC_Star, AT_MR_Star, beta, CTL );

        BLAS::Internal::SyrkUNUpdate
        ( alpha, AB_MC_Star, AB_MR_Star, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::BLAS::Internal::SyrkUNUpdateKernel
( const T alpha, const DistMatrix<T,MC,Star>& A_MC_Star,
                 const DistMatrix<T,MR,Star>& A_MR_Star,
  const T beta,        DistMatrix<T,MC,MR  >& C         )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::SyrkUNUpdateKernel");
#endif
    const Grid& grid = C.GetGrid();
#ifndef RELEASE
    if( A_MC_Star.GetGrid() != A_MR_Star.GetGrid() || 
        A_MR_Star.GetGrid() != C.GetGrid()           )
    {
        throw "A and C must be distributed over the same grid.";
    }
    if( A_MC_Star.Height() != C.Height() ||
        A_MR_Star.Height() != C.Width()  ||
        A_MC_Star.Height() != A_MR_Star.Height() ||
        A_MC_Star.Width()  != A_MR_Star.Width()    )
    {
        ostringstream msg;
        msg << "Nonconformal SyrkUNUpdateKernel: " << endl
            << "  A[MC,* ] ~ " << A_MC_Star.Height() << " x "
                               << A_MC_Star.Width()  << endl
            << "  A[MR,* ] ~ " << A_MR_Star.Height() << " x "
                               << A_MR_Star.Width()  << endl
            << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        throw msg.str();
    }
    if( A_MC_Star.ColAlignment() != C.ColAlignment() ||
        A_MR_Star.ColAlignment() != C.RowAlignment()   )
    {
        ostringstream msg;
        msg << "Misaligned SyrkUNUpdateKernel: " << endl
            << "  A[MC,* ] ~ " << A_MC_Star.ColAlignment() << endl
            << "  A[MR,* ] ~ " << A_MR_Star.ColAlignment() << endl
            << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                  C.RowAlignment() << endl;
        throw msg.str();
    }
#endif
    DistMatrix<T,MC,Star> AT_MC_Star(grid),
                          AB_MC_Star(grid);

    DistMatrix<T,MR,Star> AT_MR_Star(grid),
                          AB_MR_Star(grid);

    DistMatrix<T,MC,MR>
        CTL(grid), CTR(grid),
        CBL(grid), CBR(grid);

    DistMatrix<T,MC,MR> DTL(grid), DBR(grid);

    const unsigned half = C.Height()/2;

    BLAS::Scal( beta, C );

    LockedPartitionDown( A_MC_Star, AT_MC_Star,
                                    AB_MC_Star, half );

    LockedPartitionDown( A_MR_Star, AT_MR_Star,
                                    AB_MR_Star, half );

    PartitionDownDiagonal( C, CTL, CTR,
                              CBL, CBR, half );

    DTL.AlignWith( CTL );
    DBR.AlignWith( CBR );
    DTL.ResizeTo( CTL.Height(), CTL.Width() );
    DBR.ResizeTo( CBR.Height(), CBR.Width() );
    //------------------------------------------------------------------------//
    BLAS::Gemm( Normal, Transpose,
                alpha, AT_MC_Star.LockedLocalMatrix(),
                       AB_MR_Star.LockedLocalMatrix(),
                (T)1,  CTR.LocalMatrix()              );

    BLAS::Gemm( Normal, Transpose,
                alpha, AT_MC_Star.LockedLocalMatrix(),
                       AT_MR_Star.LockedLocalMatrix(),
                (T)0,  DTL.LocalMatrix()              );
    DTL.MakeTrapezoidal( Left, Upper );
    BLAS::Axpy( (T)1, DTL, CTL );

    BLAS::Gemm( Normal, Transpose,
                alpha, AB_MC_Star.LockedLocalMatrix(),
                       AB_MR_Star.LockedLocalMatrix(),
                (T)0,  DBR.LocalMatrix()              );
    DBR.MakeTrapezoidal( Left, Upper );
    BLAS::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//

#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Internal::SyrkUN
( const float alpha, const DistMatrix<float,MC,MR>& A,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::SyrkUN
( const double alpha, const DistMatrix<double,MC,MR>& A,
  const double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Internal::SyrkUN
( const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::SyrkUN
( const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

