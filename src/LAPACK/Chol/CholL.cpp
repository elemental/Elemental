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
#include "ElementalLAPACKInternal.h"
using namespace std;
using namespace Elemental;

// The mainline Cholesky wraps the variant 2 algorithm
template<typename T>
void
Elemental::LAPACK::Internal::CholL
( DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::CholL");
#endif
    LAPACK::Internal::CholLVar2( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

/*
   Parallelization of Variant 2 Lower Cholesky factorization. 

   Original serial update:
   ------------------------
   A11 := A11 - A10 A10^H
   A11 := Chol(A11)
   A21 := A21 - A20 A10^H
   A21 := A21 tril(A11)^-H
   ------------------------

   Our parallel update:
   -----------------------------------------------------
   A10[* ,MR] <- A10[MC,MR]
   X11[MC,* ] := A10[MC,MR] (A01[* ,MR])^H
   A11[MC,MR] := A11[MC,MR] - (SumRow(X11[MC,* ]))[* ,MR]

   A11[* ,* ] <- A11[MC,MR]   
   A11[* ,* ] := Chol(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]

   X21[MC,* ] := A20[MC,MR] (A10[* ,MR])^H
   A21[MC,MR] := A21[MC,MR] - (SumRow(X21[MC,* ]))[* ,MR]

   A21[VC,* ] <- A21[MC,MR]
   A21[VC,* ] := A21[VC,* ] tril(A11[* ,* ])^-H
   A21[MC,MR] <- A21[VC,* ]
   -----------------------------------------------------
*/
template<typename T>
void
Elemental::LAPACK::Internal::CholLVar2
( DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::CholLVar2");
    if( A.Height() != A.Width() )
    {
        if( A.GetGrid().VCRank() == 0 )
        {
            cerr << "Can only compute Cholesky factor of square matrices."
                 << endl;
        }
        DumpCallStack();
        throw exception();
    }
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        ATL(grid), ATR(grid),   A00(grid), A01(grid), A02(grid),
        ABL(grid), ABR(grid),   A10(grid), A11(grid), A12(grid),
                                A20(grid), A21(grid), A22(grid);

    // Temporary distributions
    DistMatrix<T,Star,MR  > A10_Star_MR(grid);
    DistMatrix<T,Star,Star> A11_Star_Star(grid);
    DistMatrix<T,VC,  Star> A21_VC_Star(grid);
    DistMatrix<T,MC,  Star> X11_MC_Star(grid);
    DistMatrix<T,MC,  Star> X21_MC_Star(grid);

    // Start the algorithm
    PartitionDownDiagonal( A, ATL, ATR,
                              ABL, ABR );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal( ATL, /**/ ATR,  A00, /**/ A01, A02,
                                /*************/ /******************/
                                      /**/       A10, /**/ A11, A12,
                                 ABL, /**/ ABR,  A20, /**/ A21, A22 );

        A10_Star_MR.ConformWith( A10 );
        X11_MC_Star.AlignWith( A10 );
        X11_MC_Star.ResizeTo( A11.Height(), A11.Width() );
        X21_MC_Star.AlignWith( A20 );
        X21_MC_Star.ResizeTo( A21.Height(), A21.Width() );
        //--------------------------------------------------------------------//
        A10_Star_MR = A10;
        BLAS::Gemm( Normal, ConjugateTranspose, 
                    (T)1, A10.LockedLocalMatrix(),
                          A10_Star_MR.LockedLocalMatrix(),
                    (T)0, X11_MC_Star.LocalMatrix()       );
        A11.ReduceScatterUpdate( (T)-1, X11_MC_Star );

        A11_Star_Star = A11;
        LAPACK::Chol( Lower, A11_Star_Star.LocalMatrix() );
        A11 = A11_Star_Star;

        BLAS::Gemm( Normal, ConjugateTranspose, 
                    (T)1, A20.LockedLocalMatrix(), 
                          A10_Star_MR.LockedLocalMatrix(),
                    (T)0, X21_MC_Star.LocalMatrix()       );
        A21.ReduceScatterUpdate( (T)-1, X21_MC_Star );

        A21_VC_Star = A21;
        BLAS::Trsm( Right, Lower, ConjugateTranspose, NonUnit,
                    (T)1, A11_Star_Star.LockedLocalMatrix(),
                          A21_VC_Star.LocalMatrix()           );
        A21 = A21_VC_Star;
        //--------------------------------------------------------------------//
        A10_Star_MR.FreeConstraints();
        X11_MC_Star.FreeConstraints();
        X21_MC_Star.FreeConstraints();

        SlidePartitionDownDiagonal( ATL, /**/ ATR,  A00, A01, /**/ A02,
                                         /**/       A10, A11, /**/ A12,
                                   /*************/ /******************/
                                    ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

/*
   Parallelization of Variant 3 Lower Cholesky factorization. 

   Original serial update:
   ------------------------
   A11 := Chol(A11) 
   A21 := A21 tril(A11)^-H
   A22 := A22 - A21 A21^H
   ------------------------

   Corresponding parallel update:
   -----------------------------------------------------
   A11[* ,* ] <- A11[MC,MR] 
   A11[* ,* ] := Chol(A11[* ,* ])
   A11[MC,MR] <- A11[* ,* ]
   
   A21[VC,* ] <- A21[MC,MR]
   A21[VC,* ] := A21[VC,* ] tril(A11[* ,* ])^-H
   
   A21[VR,* ] <- A21[VC,* ]
   A21[MC,* ] <- A21[VC,* ]
   A21[MR,* ] <- A21[VR,* ]
   A22[MC,MR] := A22[MC,MR] - A21[MC,* ] (A21[MR,* ])^H
   A21[MC,MR] <- A21[MC,* ]
   -----------------------------------------------------
*/
template<typename T>
void
Elemental::LAPACK::Internal::CholLVar3
( DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::CholLVar3");
    if( A.Height() != A.Width() )
    {
        if( A.GetGrid().VCRank() == 0 )
        {
            cerr << "Can only compute Cholesky factor of square matrices." 
                 << endl;
        }
        DumpCallStack();
        throw exception();
    }
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        ATL(grid), ATR(grid),  A00(grid), A01(grid), A02(grid),
        ABL(grid), ABR(grid),  A10(grid), A11(grid), A12(grid),
                               A20(grid), A21(grid), A22(grid);

    // Temporary matrix distributions
    DistMatrix<T,Star,Star> A11_Star_Star(grid);
    DistMatrix<T,VC,  Star> A21_VC_Star(grid);
    DistMatrix<T,VR,  Star> A21_VR_Star(grid);
    DistMatrix<T,MC,  Star> A21_MC_Star(grid);
    DistMatrix<T,MR,  Star> A21_MR_Star(grid);

    // Start the algorithm
    PartitionDownDiagonal( A, ATL, ATR,
                              ABL, ABR );
    while( ABR.Height() > 0 )
    {
        RepartitionDownDiagonal( ATL, /**/ ATR,  A00, /**/ A01, A02,
                                /*************/ /******************/   
                                      /**/       A10, /**/ A11, A12,
                                 ABL, /**/ ABR,  A20, /**/ A21, A22 );

        A21_MC_Star.AlignWith( A22 );
        A21_MR_Star.AlignWith( A22 );
        A21_VC_Star.AlignWith( A22 );
        A21_VR_Star.AlignWith( A22 );
        //--------------------------------------------------------------------//
        A11_Star_Star = A11;
        LAPACK::Chol( Lower, A11_Star_Star.LocalMatrix() );
        A11 = A11_Star_Star;

        A21_VC_Star = A21;
        BLAS::Trsm( Right, Lower, ConjugateTranspose, NonUnit, 
                    (T)1, A11_Star_Star.LockedLocalMatrix(), 
                          A21_VC_Star.LocalMatrix()           );

        A21_VR_Star = A21_VC_Star;
        A21_MC_Star = A21_VC_Star;
        A21_MR_Star = A21_VR_Star;
        BLAS::Internal::HerkLNUpdate
        ( (T)-1, A21_MC_Star, A21_MR_Star,(T)1, A22 );
        A21 = A21_MC_Star;
        //--------------------------------------------------------------------//
        A21_MC_Star.FreeConstraints();
        A21_MR_Star.FreeConstraints();
        A21_VC_Star.FreeConstraints();
        A21_VR_Star.FreeConstraints();

        SlidePartitionDownDiagonal( ATL, /**/ ATR,  A00, A01, /**/ A02,
                                         /**/       A10, A11, /**/ A12,
                                   /*************/ /******************/
                                    ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
} 

template void Elemental::LAPACK::Internal::CholL
( DistMatrix<float,MC,MR>& A );

template void Elemental::LAPACK::Internal::CholLVar2
( DistMatrix<float,MC,MR>& A );

template void Elemental::LAPACK::Internal::CholLVar3
( DistMatrix<float,MC,MR>& A );

template void Elemental::LAPACK::Internal::CholL
( DistMatrix<double,MC,MR>& A );

template void Elemental::LAPACK::Internal::CholLVar2
( DistMatrix<double,MC,MR>& A );

template void Elemental::LAPACK::Internal::CholLVar3
( DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void Elemental::LAPACK::Internal::CholL
( DistMatrix<scomplex,MC,MR>& A );

template void Elemental::LAPACK::Internal::CholLVar2
( DistMatrix<scomplex,MC,MR>& A );

template void Elemental::LAPACK::Internal::CholLVar3
( DistMatrix<scomplex,MC,MR>& A );

template void Elemental::LAPACK::Internal::CholL
( DistMatrix<dcomplex,MC,MR>& A );

template void Elemental::LAPACK::Internal::CholLVar2
( DistMatrix<dcomplex,MC,MR>& A );

template void Elemental::LAPACK::Internal::CholLVar3
( DistMatrix<dcomplex,MC,MR>& A );
#endif

