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

// Right Lower (Conjugate)Transpose (Non)Unit Trsm
//   X := X tril(L)^-T,
//   X := X tril(L)^-H,
//   X := X trilu(L)^-T, or
//   X := X trilu(L)^-H
template<typename T>
void
BLAS::Internal::TrsmRLT
( const Orientation orientation, 
  const Diagonal diagonal,
  const T alpha, 
  const DistMatrix<T,MC,MR>& L,
        DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::TrsmRLT");
    if( L.GetGrid() != X.GetGrid() )
        throw "L and X must be distributed over the same grid.";
    if( orientation == Normal )
        throw "TrsmRLT expects a (Conjugate)Transpose option.";
    if( L.Height() != L.Width() || X.Width() != L.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrsmRLT: " << endl
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
    DistMatrix<T,Star,Star> L11_Star_Star(grid);
    DistMatrix<T,MR,  Star> L21_MR_Star(grid);
    DistMatrix<T,MC,  Star> X1_MC_Star(grid);
    DistMatrix<T,VC,  Star> X1_VC_Star(grid);

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

        RepartitionRight( XL, /**/     XR,
                          X0, /**/ X1, X2 );

        X1_MC_Star.AlignWith( X2 );
        L21_MR_Star.AlignWith( X2 );
        //--------------------------------------------------------------------//
        L11_Star_Star = L11; // L11[*,*] <- L11[MC,MR]
        X1_VC_Star    = X1;  // X1[VC,*] <- X1[MC,MR]
        
        // X1[VC,*] := X1[VC,*] (L11[*,*])^-(T/H)
        BLAS::Trsm( Right, Lower, orientation, diagonal,
                    (T)1, L11_Star_Star.LockedLocalMatrix(),
                          X1_VC_Star.LocalMatrix()          );

        X1_MC_Star  = X1_VC_Star; // X1[MC,*]  <- X1[VC,*]
        X1          = X1_MC_Star; // X1[MC,MR] <- X1[MC,*]
        L21_MR_Star = L21;        // L21[MR,*] <- L21[MC,MR]

        // X2[MC,MR] -= X1[MC,*] (L21[MR,*])^(T/H)
        //            = X1[MC,*] (L21^(T/H))[*,MR]
        BLAS::Gemm( Normal, orientation, 
                    (T)-1, X1_MC_Star.LockedLocalMatrix(),
                           L21_MR_Star.LockedLocalMatrix(),
                    (T) 1, X2.LocalMatrix()                );
        //--------------------------------------------------------------------//
        X1_MC_Star.FreeConstraints();
        L21_MR_Star.FreeConstraints();

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

template void Elemental::BLAS::Internal::TrsmRLT
( const Orientation orientation, 
  const Diagonal diagonal,
  const float alpha, 
  const DistMatrix<float,MC,MR>& L,
        DistMatrix<float,MC,MR>& X );

template void Elemental::BLAS::Internal::TrsmRLT
( const Orientation orientation, 
  const Diagonal diagonal,
  const double alpha, 
  const DistMatrix<double,MC,MR>& L,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Internal::TrsmRLT
( const Orientation orientation, 
  const Diagonal diagonal,
  const scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& L,
        DistMatrix<scomplex,MC,MR>& X );

template void Elemental::BLAS::Internal::TrsmRLT
( const Orientation orientation, 
  const Diagonal diagonal,
  const dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& L,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

