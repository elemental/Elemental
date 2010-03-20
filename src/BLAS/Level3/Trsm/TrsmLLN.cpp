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

// Left Lower Normal (Non)Unit Trsm 
//   X := tril(L)^-1  X, or
//   X := trilu(L)^-1 X
template<typename T>
void
Elemental::BLAS::Internal::TrsmLLN
( const Diagonal diagonal,
  const T alpha, const DistMatrix<T,MC,MR>& L,
                       DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::TrsmLLN");
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
    if( L.Height() != L.Width() || L.Width() != X.Height() )
    {
        if( grid.VCRank() == 0 )
        {
            cerr << "Nonconformal TrsmLLN: " <<
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

    DistMatrix<T,MC,MR> XT(grid),  X0(grid), 
                        XB(grid),  X1(grid),
                                   X2(grid);

    // Temporary distributions
    DistMatrix<T,Star,Star> L11_Star_Star(grid);
    DistMatrix<T,MC,  Star> L21_MC_Star(grid);
    DistMatrix<T,Star,MR  > X1_Star_MR(grid);
    DistMatrix<T,Star,VR  > X1_Star_VR(grid);

    // Start the algorithm
    BLAS::Scal( alpha, X );
    LockedPartitionDownDiagonal( L, LTL, LTR,
                                    LBL, LBR );
    PartitionDown( X, XT,
                      XB );
    while( XB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal( LTL, /**/ LTR,  L00, /**/ L01, L02,
                                      /*************/ /******************/
                                            /**/       L10, /**/ L11, L12,
                                       LBL, /**/ LBR,  L20, /**/ L21, L22 );

        RepartitionDown( XT,  X0,
                        /**/ /**/
                              X1,
                         XB,  X2 );

        L21_MC_Star.AlignWith( X2 );
        X1_Star_MR.AlignWith( X2 );
        //--------------------------------------------------------------------//
        L11_Star_Star = L11; // L11[*,*] <- L11[MC,MR]
        X1_Star_VR    = X1;  // X1[*,VR] <- X1[MC,MR]

        // X1[*,VR] := (L11[*,*])^-1 X1[*,VR]
        BLAS::Trsm( Left, Lower, Normal, diagonal,
                    (T)1, L11_Star_Star.LockedLocalMatrix(),
                          X1_Star_VR.LocalMatrix()          );

        X1_Star_MR  = X1_Star_VR; // X1[*,MR]  <- X1[*,VR]
        X1          = X1_Star_MR; // X1[MC,MR] <- X1[*,MR]
        L21_MC_Star = L21;        // L21[MC,*] <- L21[MC,MR]
        
        // X2[MC,MR] -= L21[MC,*] X1[*,MR]
        BLAS::Gemm( Normal, Normal, 
                    (T)-1, L21_MC_Star.LockedLocalMatrix(),
                           X1_Star_MR.LockedLocalMatrix(),
                    (T) 1, X2.LocalMatrix()                );
        //--------------------------------------------------------------------//
        L21_MC_Star.FreeConstraints();
        X1_Star_MR.FreeConstraints();

        SlideLockedPartitionDownDiagonal( LTL, /**/ LTR,  L00, L01, /**/ L02,
                                               /**/       L10, L11, /**/ L12, 
                                         /*************/ /******************/
                                          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionDown( XT,  X0,
                                 X1,
                           /**/ /**/
                            XB,  X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Internal::TrsmLLN
( const Diagonal diagonal,
  const float alpha, 
  const DistMatrix<float,MC,MR>& L,
        DistMatrix<float,MC,MR>& X );

template void Elemental::BLAS::Internal::TrsmLLN
( const Diagonal diagonal,
  const double alpha, 
  const DistMatrix<double,MC,MR>& L,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Internal::TrsmLLN
( const Diagonal diagonal,
  const scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& L,
        DistMatrix<scomplex,MC,MR>& X );

template void Elemental::BLAS::Internal::TrsmLLN
( const Diagonal diagonal,
  const dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& L,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

