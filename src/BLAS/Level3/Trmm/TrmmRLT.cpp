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

// Right Lower (Conjugate)Transpose (Non)Unit Trmm
//   X := X tril(L)^T,
//   X := X tril(L)^H,
//   X := X trilu(L)^T, or
//   X := X trilu(L)^H
template<typename T>
void
BLAS::Internal::TrmmRLT
( const Orientation orientation, 
  const Diagonal diagonal,
  const T alpha, const DistMatrix<T,MC,MR>& L,
                       DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::TrmmRLT");
    if( L.GetGrid() != X.GetGrid() )
        throw "L and X must be distributed over the same grid.";
    if( orientation == Normal )
        throw "TrmmRLT expects a (Conjugate)Transpose option.";
    if( L.Height() != L.Width() || X.Width() != L.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrmmRLT: " << endl
            << "  L ~ " << L.Height() << " x " << L.Width() << endl
            << "  X ~ " << X.Height() << " x " << X.Width() << endl;
        throw msg.str();
    }
#endif
    const Grid& grid = L.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        LTL(grid), LTR(grid),  L00(grid), L01(grid), L02(grid),
        LBL(grid), LBR(grid),  L10(grid), L11(grid), L12(grid),
                               L20(grid), L21(grid), L22(grid);

    DistMatrix<T,MC,MR> XL(grid), XR(grid),
                        X0(grid), X1(grid), X2(grid);

    // Temporary distributions
    DistMatrix<T,Star,MR  > L10_Star_MR(grid);
    DistMatrix<T,Star,Star> L11_Star_Star(grid);
    DistMatrix<T,VC,  Star> X1_VC_Star(grid);
    DistMatrix<T,MC,  Star> D1_MC_Star(grid);

    // Start the algorithm
    BLAS::Scal( alpha, X );
    LockedPartitionUpDiagonal( L, LTL, LTR,
                                  LBL, LBR );
    PartitionLeft( X, XL, XR );
    while( XL.Width() > 0 )
    {
        LockedRepartitionUpDiagonal( LTL, /**/ LTR,  L00, L01, /**/ L02,
                                          /**/       L10, L11, /**/ L12,
                                    /*************/ /******************/
                                     LBL, /**/ LBR,  L20, L21, /**/ L22 );

        RepartitionLeft( XL,     /**/ XR,
                         X0, X1, /**/ X2 );

        L10_Star_MR.ConformWith( X0 );
        D1_MC_Star.AlignWith( X1 );
        D1_MC_Star.ResizeTo( X1.Height(), X1.Width() );
        //--------------------------------------------------------------------//
        X1_VC_Star = X1;
        L11_Star_Star = L11;
        BLAS::Trmm( Right, Lower, orientation, diagonal,
                    (T)1, L11_Star_Star.LockedLocalMatrix(),
                          X1_VC_Star.LocalMatrix()          );
        X1 = X1_VC_Star;
 
        L10_Star_MR = L10;
        BLAS::Gemm( Normal, orientation, 
                    (T)1, X0.LockedLocalMatrix(),
                          L10_Star_MR.LockedLocalMatrix(),
                    (T)0, D1_MC_Star.LocalMatrix()        );
        X1.ReduceScatterUpdate( (T)1, D1_MC_Star );
       //--------------------------------------------------------------------//
        L10_Star_MR.FreeConstraints();
        D1_MC_Star.FreeConstraints();

        SlideLockedPartitionUpDiagonal( LTL, /**/ LTR,  L00, /**/ L01, L02,
                                       /*************/ /******************/
                                             /**/       L10, /**/ L11, L12,
                                        LBL, /**/ LBR,  L20, /**/ L21, L22 );

        SlidePartitionLeft( XL, /**/     XR,
                            X0, /**/ X1, X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Internal::TrmmRLT
( const Orientation orientation, 
  const Diagonal diagonal,
  const float alpha, 
  const DistMatrix<float,MC,MR>& L,
        DistMatrix<float,MC,MR>& X );

template void Elemental::BLAS::Internal::TrmmRLT
( const Orientation orientation, 
  const Diagonal diagonal,
  const double alpha, 
  const DistMatrix<double,MC,MR>& L,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Internal::TrmmRLT
( const Orientation orientation, 
  const Diagonal diagonal,
  const scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& L,
        DistMatrix<scomplex,MC,MR>& X );

template void Elemental::BLAS::Internal::TrmmRLT
( const Orientation orientation, 
  const Diagonal diagonal,
  const dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& L,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

