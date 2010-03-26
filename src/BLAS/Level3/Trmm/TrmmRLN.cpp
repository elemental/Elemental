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

// Right Lower Normal (Non)Unit Trmm
//   X := X tril(L), and
//   X := X trilu(L)
template<typename T>
void
BLAS::Internal::TrmmRLN
( const Diagonal diagonal,
  const T alpha, const DistMatrix<T,MC,MR>& L,
                       DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::TrmmRLN");
#endif
    const Grid& grid = L.GetGrid();
#ifndef RELEASE
    if( L.GetGrid() != X.GetGrid() )
    {
        if( grid.VCRank() == 0 )
            cerr << "L and X must be distributed over the same grid." << endl;
        DumpCallStack();
        throw exception();
    }
    if( L.Height() != L.Width() || X.Width() != L.Height() )
    {
        if( grid.VCRank() == 0 )
        {
            cerr << "Nonconformal TrmmRLN: " <<
            endl << "  L ~ " << L.Height() << " x " << L.Width() <<
            endl << "  X ~ " << X.Height() << " x " << X.Width() << endl;
        }
        DumpCallStack();
        throw exception();
    }
#endif
    // Matrix views
    DistMatrix<T,MC,MR> 
        LTL(grid), LTR(grid),  L00(grid), L01(grid), L02(grid),
        LBL(grid), LBR(grid),  L10(grid), L11(grid), L12(grid),
                               L20(grid), L21(grid), L22(grid);

    DistMatrix<T,MC,MR> XL(grid), XR(grid),
                        X0(grid), X1(grid), X2(grid);

    // Temporary distributions
    DistMatrix<T,Star,Star> L11_Star_Star(grid);
    DistMatrix<T,MR,  Star> L21_MR_Star(grid);
    DistMatrix<T,VC,  Star> X1_VC_Star(grid);
    DistMatrix<T,MC,  Star> D1_MC_Star(grid);

    // Start the algorithm
    BLAS::Scal( alpha, X );
    LockedPartitionDownDiagonal( L, LTL, LTR,
                                    LBL, LBR );
    PartitionRight( X, XL, XR );
    while( XR.Width() > 0 )
    {
        LockedRepartitionDownDiagonal( LTL, /**/ LTR,  L00, /**/ L01, L02,
                                      /*************/ /******************/
                                            /**/       L10, /**/ L11, L12,
                                       LBL, /**/ LBR,  L20, /**/ L21, L22 );
 
        RepartitionRight( XL, /**/ XR,
                          X0, /**/ X1, X2 );

        L21_MR_Star.ConformWith( X2 );
        D1_MC_Star.AlignWith( X1 );
        D1_MC_Star.ResizeTo( X1.Height(), X1.Width() );
        //--------------------------------------------------------------------//
        X1_VC_Star = X1;
        L11_Star_Star = L11;
        BLAS::Trmm( Right, Lower, Normal, diagonal,
                    (T)1, L11_Star_Star.LockedLocalMatrix(),
                          X1_VC_Star.LocalMatrix()          );
        X1 = X1_VC_Star;
 
        L21_MR_Star = L21;
        BLAS::Gemm( Normal, Normal, 
                    (T)1, X2.LockedLocalMatrix(),
                          L21_MR_Star.LockedLocalMatrix(),
                    (T)0, D1_MC_Star.LocalMatrix()        );
        X1.ReduceScatterUpdate( (T)1, D1_MC_Star );
       //--------------------------------------------------------------------//
        L21_MR_Star.FreeConstraints();
        D1_MC_Star.FreeConstraints();

        SlideLockedPartitionDownDiagonal( LTL, /**/ LTR,  L00, L01, /**/ L02,
                                               /**/       L10, L11, /**/ L12,
                                         /*************/ /******************/
                                          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionRight( XL,     /**/ XR,
                             X0, X1, /**/ X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Internal::TrmmRLN
( const Diagonal diagonal,
  const float alpha, const DistMatrix<float,MC,MR>& L,
                           DistMatrix<float,MC,MR>& X );

template void Elemental::BLAS::Internal::TrmmRLN
( const Diagonal diagonal,
  const double alpha, const DistMatrix<double,MC,MR>& L,
                            DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Internal::TrmmRLN
( const Diagonal diagonal,
  const scomplex alpha, const DistMatrix<scomplex,MC,MR>& L,
                              DistMatrix<scomplex,MC,MR>& X );

template void Elemental::BLAS::Internal::TrmmRLN
( const Diagonal diagonal,
  const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& L,
                              DistMatrix<dcomplex,MC,MR>& X );
#endif

