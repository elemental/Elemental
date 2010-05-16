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
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::blas::internal::GemmTN
( Orientation orientationOfA,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmTN");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "{A,B,C} must be distributed over the same grid.";
    if( orientationOfA == Normal )
        throw "GemmTN assumes A is (Conjugate)Transposed.";
#endif
    const int m = C.Height();
    const int n = C.Width();
    const int k = A.Height();
    const float weightTowardsC = 2.0;

    if( m <= n && weightTowardsC*m <= k )
    {
#ifndef RELEASE
        if( A.GetGrid().VCRank() == 0 )
            cout << "  GemmTN routing to GemmTNB." << endl;
#endif
        blas::internal::GemmTNB( orientationOfA, alpha, A, B, beta, C );
    }
    else if( n <= m && weightTowardsC*n <= k )
    {
#ifndef RELEASE
        if( A.GetGrid().VCRank() == 0 )
            cout << "  GemmTN routing to GemmTNA." << endl;
#endif
        blas::internal::GemmTNA( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
#ifndef RELEASE
        if( A.GetGrid().VCRank() == 0 )
            cout << "  GemmTN routing to GemmTNC." << endl;
#endif
        blas::internal::GemmTNC( orientationOfA, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Transpose Normal Gemm that avoids communicating the matrix A.
template<typename T> 
void
elemental::blas::internal::GemmTNA
( Orientation orientationOfA,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmTNA");    
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "{A,B,C} must be distributed over the same grid.";
    if( orientationOfA == Normal )
        throw "GemmTNA assumes A is (Conjugate)Transposed.";
    if( A.Width()  != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Height() != B.Height()   )
    {
        ostringstream msg;
        msg << "Nonconformal GemmTNA: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B ~ " << B.Height() << " x " << B.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        const string& s = msg.str();
        throw s.c_str();
    }
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> BL(grid), BR(grid),
                        B0(grid), B1(grid), B2(grid);

    DistMatrix<T,MC,MR> CL(grid), CR(grid),
                        C0(grid), C1(grid), C2(grid);

    // Temporary distributions
    DistMatrix<T,MC,Star> B1_MC_Star(grid);
    DistMatrix<T,MR,Star> D1_MR_Star(grid);
    DistMatrix<T,MR,MC  > D1_MR_MC(grid);
    DistMatrix<T,MC,MR  > D1(grid);

    // Start the algorithm
    blas::Scal( beta, C );
    LockedPartitionRight( B, BL, BR );
    PartitionRight( C, CL, CR );
    while( BR.Width() > 0 )
    {
        LockedRepartitionRight( BL, /**/     BR,
                                B0, /**/ B1, B2 );
 
        RepartitionRight( CL, /**/     CR,
                          C0, /**/ C1, C2 );

        B1_MC_Star.ConformWith( A );
        D1_MR_Star.AlignWith( A );
        D1_MR_Star.ResizeTo( C1.Height(), C1.Width() );
        D1.AlignWith( C1 );
        //--------------------------------------------------------------------//
        B1_MC_Star = B1; // B1[MC,*] <- B1[MC,MR]

        // D1[MR,*] := alpha (A1[MC,MR])^T B1[MC,*]
        //           = alpha (A1^T)[MR,MC] B1[MC,*]
        blas::Gemm( orientationOfA, Normal, 
                    alpha, A.LockedLocalMatrix(),
                           B1_MC_Star.LockedLocalMatrix(),
                    (T)0,  D1_MR_Star.LocalMatrix()       );

        // C1[MC,MR] += scattered & transposed D1[MR,*] summed over grid cols
        D1_MR_MC.ReduceScatterFrom( D1_MR_Star );
        D1 = D1_MR_MC; 
        blas::Axpy( (T)1, D1, C1 );
        //--------------------------------------------------------------------//
        B1_MC_Star.FreeConstraints();
        D1_MR_Star.FreeConstraints();
        D1.FreeConstraints();

        SlideLockedPartitionRight( BL,     /**/ BR,
                                   B0, B1, /**/ B2 );

        SlidePartitionRight( CL,     /**/ CR,
                             C0, C1, /**/ C2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Transpose Normal Gemm that avoids communicating the matrix B.
template<typename T> 
void
elemental::blas::internal::GemmTNB
( Orientation orientationOfA,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmTNB");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "{A,B,C} must be distributed over the same grid.";
    if( orientationOfA == Normal )
        throw "GemmTNB assumes A is (Conjugate)Transposed.";
    if( A.Width()  != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Height() != B.Height()   )
    {
        ostringstream msg;
        msg << "Nonconformal GemmTNB: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B ~ " << B.Height() << " x " << B.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        const string& s = msg.str();
        throw s.c_str();
    }
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> AL(grid), AR(grid),
                        A0(grid), A1(grid), A2(grid);

    DistMatrix<T,MC,MR> CT(grid),  C0(grid),
                        CB(grid),  C1(grid),
                                   C2(grid);

    // Temporary distributions
    DistMatrix<T,MC,Star> A1_MC_Star(grid);
    DistMatrix<T,Star,MR> D1_Star_MR(grid);

    // Start the algorithm
    blas::Scal( beta, C );
    LockedPartitionRight( A, AL, AR );
    PartitionDown( C, CT,
                      CB );
    while( AR.Width() > 0 )
    {
        LockedRepartitionRight( AL, /**/     AR,
                                A0, /**/ A1, A2 );

        RepartitionDown( CT,  C0,
                        /**/ /**/
                              C1,
                         CB,  C2 );

        A1_MC_Star.ConformWith( B );
        D1_Star_MR.AlignWith( B );
        D1_Star_MR.ResizeTo( C1.Height(), C1.Width() );
        //--------------------------------------------------------------------//
        A1_MC_Star = A1; // A1[MC,*] <- A1[MC,MR]

        // D1[*,MR] := alpha (A1[MC,*])^T B[MC,MR]
        //           = alpha (A1^T)[*,MC] B[MC,MR]
        blas::Gemm( orientationOfA, Normal, 
                    alpha, A1_MC_Star.LockedLocalMatrix(),
                           B.LockedLocalMatrix(),
                    (T)0,  D1_Star_MR.LocalMatrix()       );

        // C1[MC,MR] += scattered result of D1[*,MR] summed over grid cols
        C1.ReduceScatterUpdate( (T)1, D1_Star_MR );
        //--------------------------------------------------------------------//
        A1_MC_Star.FreeConstraints();
        D1_Star_MR.FreeConstraints();

        SlideLockedPartitionRight( AL,     /**/ AR,
                                   A0, A1, /**/ A2 );

        SlidePartitionDown( CT,  C0,
                                 C1,
                           /**/ /**/
                            CB,  C2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Transpose Normal Gemm that avoids communicating the matrix C.
template<typename T> 
void
elemental::blas::internal::GemmTNC
( Orientation orientationOfA,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmTNC");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "{A,B,C} must be distributed over the same grid.";
    if( orientationOfA == Normal )
        throw "GemmTNC assumes A is (Conjugate)Transposed.";
    if( A.Width()  != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Height() != B.Height()   )
    {
        ostringstream msg;
        msg << "Nonconformal GemmTNC: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B ~ " << B.Height() << " x " << B.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        const string& s = msg.str();
        throw s.c_str();
    }
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> AT(grid),  A0(grid),
                        AB(grid),  A1(grid),
                                   A2(grid);

    DistMatrix<T,MC,MR> BT(grid),  B0(grid),
                        BB(grid),  B1(grid),
                                   B2(grid);

    // Temporary distributions
    DistMatrix<T,Star,MC> A1_Star_MC(grid);
    DistMatrix<T,Star,MR> B1_Star_MR(grid);

    // Start the algorithm
    blas::Scal( beta, C );
    LockedPartitionDown( A, AT,
                            AB );
    LockedPartitionDown( B, BT,
                            BB );
    while( AB.Height() > 0 )
    {
        LockedRepartitionDown( AT,  A0,
                              /**/ /**/
                                    A1,
                               AB,  A2 );

        LockedRepartitionDown( BT,  B0,
                              /**/ /**/
                                    B1,
                               BB,  B2 );

        A1_Star_MC.AlignWith( C );
        B1_Star_MR.AlignWith( C );
        //--------------------------------------------------------------------//
        A1_Star_MC = A1; // A1[*,MC] <- A1[MC,MR]
        B1_Star_MR = B1; // B1[*,MR] <- B1[MC,MR]

        // C[MC,MR] += alpha (A1[*,MC])^T B1[*,MR]
        //           = alpha (A1^T)[MC,*] B1[*,MR]
        blas::Gemm( orientationOfA, Normal, 
                    alpha, A1_Star_MC.LockedLocalMatrix(),
                           B1_Star_MR.LockedLocalMatrix(),
                    (T)1,  C.LocalMatrix()                );
        //--------------------------------------------------------------------//
        A1_Star_MC.FreeConstraints();
        B1_Star_MR.FreeConstraints();

        SlideLockedPartitionDown( AT,  A0,
                                       A1,
                                 /**/ /**/
                                  AB,  A2 );

        SlideLockedPartitionDown( BT,  B0,
                                       B1,
                                 /**/ /**/
                                  BB,  B2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::GemmTN
( Orientation orientationOfA,
  float alpha, const DistMatrix<float,MC,MR>& A,      
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::GemmTN
( Orientation orientationOfA,
  double alpha, const DistMatrix<double,MC,MR>& A,         
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::GemmTN
( Orientation orientationOfA,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A, 
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::GemmTN
( Orientation orientationOfA,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A, 
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif
