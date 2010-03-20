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

// Right Upper (Conjugate)Transpose (Non)Unit Trmm
//   X := X triu(U)^T, 
//   X := X triu(U)^H,
//   X := X triuu(U)^T, or
//   X := X triuu(U)^H
template<typename T>
void
BLAS::Internal::TrmmRUT
( const Orientation orientation, 
  const Diagonal diagonal,
  const T alpha, const DistMatrix<T,MC,MR>& U,
                       DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::TrmmRUT");
#endif
    const Grid& grid = U.GetGrid();
#ifndef RELEASE
    if( U.GetGrid() != X.GetGrid() )
    {
        if( grid.VCRank() == 0 )
            cerr << "U and X must be distributed over the same grid." << endl;
        DumpCallStack();
        throw exception();
    }
    if( orientation == Normal )
    {
        if( grid.VCRank() == 0 )
            cerr << "TrmmRUT expects a (Conjugate)Transpose option." << endl;
        DumpCallStack();
        throw exception();
    }
    if( U.Height() != U.Width() || X.Width() != U.Height() )
    {
        if( grid.VCRank() == 0 )
        {
            cerr << "Nonconformal TrmmRUT: " <<
            endl << "  U ~ " << U.Height() << " x " << U.Width() <<
            endl << "  X ~ " << X.Height() << " x " << X.Width() << endl;
        }
        DumpCallStack();
        throw exception();
    }
#endif
    // Matrix views
    DistMatrix<T,MC,MR> 
        UTL(grid), UTR(grid),  U00(grid), U01(grid), U02(grid),
        UBL(grid), UBR(grid),  U10(grid), U11(grid), U12(grid),
                               U20(grid), U21(grid), U22(grid);

    DistMatrix<T,MC,MR> XL(grid), XR(grid),
                        X0(grid), X1(grid), X2(grid);

    // Temporary distributions
    DistMatrix<T,Star,Star> U11_Star_Star(grid);
    DistMatrix<T,Star,MR  > U12_Star_MR(grid);
    DistMatrix<T,VC,  Star> X1_VC_Star(grid);
    DistMatrix<T,MC,  Star> D1_MC_Star(grid);
    
    // Start the algorithm
    BLAS::Scal( alpha, X );
    LockedPartitionDownDiagonal( U, UTL, UTR,
                                    UBL, UBR );
    PartitionRight( X, XL, XR );
    while( XR.Width() > 0 )
    {
        LockedRepartitionDownDiagonal( UTL, /**/ UTR,  U00, /**/ U01, U02,
                                      /*************/ /******************/
                                            /**/       U10, /**/ U11, U12,
                                       UBL, /**/ UBR,  U20, /**/ U21, U22 );

        RepartitionRight( XL, /**/ XR,
                          X0, /**/ X1, X2 );

        U12_Star_MR.ConformWith( X2 );
        D1_MC_Star.AlignWith( X1 );
        D1_MC_Star.ResizeTo( X1.Height(), X1.Width() );
        //--------------------------------------------------------------------//
        X1_VC_Star = X1;
        U11_Star_Star = U11;
        BLAS::Trmm( Right, Upper, orientation, diagonal,
                    (T)1, U11_Star_Star.LockedLocalMatrix(),
                          X1_VC_Star.LocalMatrix()          );
        X1 = X1_VC_Star;
 
        U12_Star_MR = U12;
        BLAS::Gemm( Normal, orientation, 
                    (T)1, X2.LockedLocalMatrix(),
                          U12_Star_MR.LockedLocalMatrix(),
                    (T)0, D1_MC_Star.LocalMatrix()        );
        X1.ReduceScatterUpdate( (T)1, D1_MC_Star );
       //--------------------------------------------------------------------//
        U12_Star_MR.FreeConstraints();
        D1_MC_Star.FreeConstraints();

        SlideLockedPartitionDownDiagonal( UTL, /**/ UTR,  U00, U01, /**/ U02,
                                               /**/       U10, U11, /**/ U12, 
                                         /*************/ /******************/
                                          UBL, /**/ UBR,  U20, U21, /**/ U22 );

        SlidePartitionRight( XL,     /**/ XR,
                             X0, X1, /**/ X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Internal::TrmmRUT
( const Orientation orientation, 
  const Diagonal diagonal,
  const float alpha, 
  const DistMatrix<float,MC,MR>& U,
        DistMatrix<float,MC,MR>& X );

template void Elemental::BLAS::Internal::TrmmRUT
( const Orientation orientation, 
  const Diagonal diagonal,
  const double alpha, 
  const DistMatrix<double,MC,MR>& U,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Internal::TrmmRUT
( const Orientation orientation, 
  const Diagonal diagonal,
  const scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& U,
        DistMatrix<scomplex,MC,MR>& X );

template void Elemental::BLAS::Internal::TrmmRUT
( const Orientation orientation, 
  const Diagonal diagonal,
  const dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& U,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

