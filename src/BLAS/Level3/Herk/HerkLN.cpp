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
Elemental::BLAS::Internal::HerkLN
( const T alpha, const DistMatrix<T,MC,MR>& A,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::HerkLN");
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
    if( A.Height() != C.Height() || A.Height() != C.Width() )
    {
        if( grid.VCRank() == 0 )
        {
            cerr << "Nonconformal HerkLN:" <<
            endl << "  A ~ " << A.Height() << " x " << A.Width() <<
            endl << "  C ~ " << C.Height() << " x " << C.Width() <<  endl;
        }
        DumpCallStack();
        throw exception();
    }
#endif
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

        BLAS::Internal::HerkLNUpdate
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
Elemental::BLAS::Internal::HerkLNUpdate
( const T alpha, const DistMatrix<T,MC,Star>& A_MC_Star,
                 const DistMatrix<T,MR,Star>& A_MR_Star,
  const T beta,        DistMatrix<T,MC,MR  >& C         )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::HerkLNUpdate");
#endif
    const Grid& grid = C.GetGrid();
#ifndef RELEASE
    if( A_MC_Star.GetGrid() != A_MR_Star.GetGrid() || 
        A_MR_Star.GetGrid() != C.GetGrid()            )
    {
        if( grid.VCRank() == 0 )
            cerr << "A and C must be distributed over the same grid." << endl;
        DumpCallStack();
        throw exception();
    }
    if( A_MC_Star.Height() != C.Height() ||
        A_MR_Star.Height() != C.Width()  ||
        A_MC_Star.Height() != A_MR_Star.Height() ||
        A_MC_Star.Width()  != A_MR_Star.Width()    )
    {
        if( grid.VCRank() == 0 )
        {
            cerr << "Nonconformal HerkLNUpdate: " <<
            endl << "  A[MC,* ] ~ " << A_MC_Star.Height() << " x "
                                    << A_MC_Star.Width()  <<
            endl << "  A[MR,* ] ~ " << A_MR_Star.Height() << " x "
                                    << A_MR_Star.Width()  <<
            endl << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        }
        DumpCallStack();
        throw exception();
    }
    if( A_MC_Star.ColAlignment() != C.ColAlignment() ||
        A_MR_Star.ColAlignment() != C.RowAlignment()   )
    {
        if( grid.VCRank() == 0 )
        {
            cerr << "Misaligned HerkLNUpdate: " <<
            endl << "  A[MC,* ] ~ " << A_MC_Star.ColAlignment() <<
            endl << "  A[MR,* ] ~ " << A_MR_Star.ColAlignment() <<
            endl << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                       C.RowAlignment() << endl;
        }
        DumpCallStack();
        throw exception();
    }
#endif
    if( C.Height() < 2*grid.Width()*Blocksize() )
    {
        BLAS::Internal::HerkLNUpdateKernel
        ( alpha, A_MC_Star, A_MR_Star, beta, C );
    }
    else
    {
        // Split C in four roughly equal pieces, perform a large gemm on CBL
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
        ( Normal, ConjugateTranspose,
          alpha, AB_MC_Star.LockedLocalMatrix(),
                 AT_MR_Star.LockedLocalMatrix(),
          beta,  CBL.LocalMatrix()              );

        // Recurse
        BLAS::Internal::HerkLNUpdate
        ( alpha, AT_MC_Star, AT_MR_Star, beta, CTL );

        BLAS::Internal::HerkLNUpdate
        ( alpha, AB_MC_Star, AB_MR_Star, beta, CBR );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::BLAS::Internal::HerkLNUpdateKernel
( const T alpha, const DistMatrix<T,MC,Star>& A_MC_Star,
                 const DistMatrix<T,MR,Star>& A_MR_Star,
  const T beta,        DistMatrix<T,MC,MR  >& C         )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::HerkLNUpdateKernel");
#endif
    const Grid& grid = C.GetGrid();
#ifndef RELEASE
    if( A_MC_Star.GetGrid() != A_MR_Star.GetGrid() || 
        A_MR_Star.GetGrid() != C.GetGrid()            )
    {
        if( grid.VCRank() == 0 )
            cerr << "A and C must be distributed over the same grid." << endl;
        DumpCallStack();
        throw exception();
    }
    if( A_MC_Star.Height() != C.Height() ||
        A_MR_Star.Height() != C.Width()  ||
        A_MC_Star.Height() != A_MR_Star.Height() ||
        A_MC_Star.Width()  != A_MR_Star.Width()    )
    {
        if( grid.VCRank() == 0 )
        {
            cerr << "Nonconformal HerkLNUpdateKernel: " <<
            endl << "  A[MC,* ] ~ " << A_MC_Star.Height() << " x "
                                    << A_MC_Star.Width()  <<
            endl << "  A[MR,* ] ~ " << A_MR_Star.Height() << " x "
                                    << A_MR_Star.Width()  <<
            endl << "  C[MC,MR] ~ " << C.Height() << " x " << C.Width() << endl;
        }
        DumpCallStack();
        throw exception();
    }
    if( A_MC_Star.ColAlignment() != C.ColAlignment() ||
        A_MR_Star.ColAlignment() != C.RowAlignment()   )
    {
        if( grid.VCRank() == 0 )
        {
            cerr << "Misaligned HerkLNUpdateKernel: " <<
            endl << "  A[MC,* ] ~ " << A_MC_Star.ColAlignment() <<
            endl << "  A[MR,* ] ~ " << A_MR_Star.ColAlignment() <<
            endl << "  C[MC,MR] ~ " << C.ColAlignment() << " , " <<
                                       C.RowAlignment() << endl;
        }
        DumpCallStack();
        throw exception();
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
    BLAS::Gemm( Normal, ConjugateTranspose,
                alpha, AB_MC_Star.LockedLocalMatrix(),
                       AT_MR_Star.LockedLocalMatrix(),
                (T)1,  CBL.LocalMatrix()              );

    BLAS::Gemm( Normal, ConjugateTranspose,
                alpha, AT_MC_Star.LockedLocalMatrix(),
                       AT_MR_Star.LockedLocalMatrix(),
                (T)0,  DTL.LocalMatrix()              );
    DTL.MakeTrapezoidal( Left, Lower );
    BLAS::Axpy( (T)1, DTL, CTL );

    BLAS::Gemm( Normal, ConjugateTranspose,
                alpha, AB_MC_Star.LockedLocalMatrix(),
                       AB_MR_Star.LockedLocalMatrix(),
                (T)0,  DBR.LocalMatrix()              );
    DBR.MakeTrapezoidal( Left, Lower );
    BLAS::Axpy( (T)1, DBR, CBR );
    //------------------------------------------------------------------------//

#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Internal::HerkLN
( const float alpha, const DistMatrix<float,MC,MR>& A,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::HerkLNUpdate
( const float alpha, const DistMatrix<float,MC,Star>& A_MC_Star,
                     const DistMatrix<float,MR,Star>& A_MR_Star,
  const float beta,        DistMatrix<float,MC,MR  >& C         );

template void Elemental::BLAS::Internal::HerkLN
( const double alpha, const DistMatrix<double,MC,MR>& A,
  const double beta,        DistMatrix<double,MC,MR>& C );

template void Elemental::BLAS::Internal::HerkLNUpdate
( const double alpha, const DistMatrix<double,MC,Star>& A_MC_Star,
                      const DistMatrix<double,MR,Star>& A_MR_Star,
  const double beta,        DistMatrix<double,MC,MR  >& C         );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Internal::HerkLN
( const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::HerkLNUpdate
( const scomplex alpha, 
  const DistMatrix<scomplex,MC,Star>& A_MC_Star,
  const DistMatrix<scomplex,MR,Star>& A_MR_Star,
  const scomplex beta,        
        DistMatrix<scomplex,MC,MR  >& C         );

template void Elemental::BLAS::Internal::HerkLN
( const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::HerkLNUpdate
( const dcomplex alpha, 
  const DistMatrix<dcomplex,MC,Star>& A_MC_Star,
  const DistMatrix<dcomplex,MR,Star>& A_MR_Star,
  const dcomplex beta,        
        DistMatrix<dcomplex,MC,MR  >& C         );
#endif

