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
Elemental::BLAS::Internal::GemmNT
( const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::GemmNT");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "{A,B,C} must be distributed over the same grid.";
    if( orientationOfB == Normal )
        throw "GemmNT requires that B be (Conjugate)Transposed.";
#endif
    const int m = C.Height();
    const int n = C.Width();
    const int k = A.Width();
    const float weightTowardsC = 2.0;

    if( m <= n && weightTowardsC*m <= k )
    {
#ifndef RELEASE
        if( A.GetGrid().VCRank() == 0 )
            cout << "  GemmNT routing to GemmNTB." << endl;
#endif
        BLAS::Internal::GemmNTB( orientationOfB, alpha, A, B, beta, C );
    }
    else if( n <= m && weightTowardsC*n <= k )
    {
#ifndef RELEASE
        if( A.GetGrid().VCRank() == 0 )
            cout << "  GemmNT routing to GemmNTA." << endl;
#endif
        BLAS::Internal::GemmNTA( orientationOfB, alpha, A, B, beta, C );
    }
    else
    {
#ifndef RELEASE
        if( A.GetGrid().VCRank() == 0 )
            cout << "  GemmNT routing to GemmNTC." << endl;
#endif
        BLAS::Internal::GemmNTC( orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Normal Transpose Gemm that avoids communicating the matrix A.
template<typename T>
void
Elemental::BLAS::Internal::GemmNTA
( const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::GemmNTA");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "{A,B,C} must be distributed over the same grid.";
    if( orientationOfB == Normal )
        throw "GemmTNA requires that B be (Conjugate)Transposed.";
    if( A.Height() != C.Height() ||
        B.Height() != C.Width()  ||
        A.Width()  != B.Width()    )
    {
        ostringstream msg;
        msg << "Nonconformal GemmNTA: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B ~ " << B.Height() << " x " << B.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        throw msg.str();
    }
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> BT(grid),  B0(grid),
                        BB(grid),  B1(grid),
                                   B2(grid);

    DistMatrix<T,MC,MR> CL(grid), CR(grid),
                        C0(grid), C1(grid), C2(grid);

    // Temporary distributions
    DistMatrix<T,Star,MR> B1_Star_MR(grid);
    DistMatrix<T,MC,Star> D1_MC_Star(grid);

    // Start the algorithm
    BLAS::Scal( beta, C );
    LockedPartitionDown( B, BT,
                            BB );
    PartitionRight( C, CL, CR );
    while( BB.Height() > 0 )
    {
        LockedRepartitionDown( BT,  B0,
                              /**/ /**/
                                    B1,
                               BB,  B2 );

        RepartitionRight( CL, /**/     CR,
                          C0, /**/ C1, C2 );

        B1_Star_MR.ConformWith( A );
        D1_MC_Star.AlignWith( A );
        D1_MC_Star.ResizeTo( C1.Height(), C1.Width() );
        //--------------------------------------------------------------------//
        B1_Star_MR = B1; // B1[*,MR] <- B1[MC,MR]

        // C1[MC,*] := alpha A[MC,MR] (B1[*,MR])^T
        //           = alpha A[MC,MR] (B1^T)[MR,*]
        BLAS::Gemm( Normal, orientationOfB, 
                    alpha, A.LockedLocalMatrix(),
                           B1_Star_MR.LockedLocalMatrix(),
                    (T)0,  D1_MC_Star.LocalMatrix()       );

        // C1[MC,MR] += scattered result of D1[MC,*] summed over grid rows
        C1.ReduceScatterUpdate( (T)1, D1_MC_Star );
        //--------------------------------------------------------------------//
        B1_Star_MR.FreeConstraints();
        D1_MC_Star.FreeConstraints();

        SlideLockedPartitionDown( BT,  B0,
                                       B1,
                                 /**/ /**/
                                  BB,  B2 );

        SlidePartitionRight( CL,     /**/ CR,
                             C0, C1, /**/ C2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Normal Transpose Gemm that avoids communicating the matrix B.
template<typename T>
void
Elemental::BLAS::Internal::GemmNTB
( const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::GemmNTB");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "{A,B,C} must be distributed over the same grid.";
    if( orientationOfB == Normal )
        throw "GemmNTB requires that B be (Conjugate)Transposed.";
    if( A.Height() != C.Height() ||
        B.Height() != C.Width()  ||
        A.Width()  != B.Width()    )
    {
        ostringstream msg;
        msg << "Nonconformal GemmNTB: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B ~ " << B.Height() << " x " << B.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        throw msg.str();
    }
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> AT(grid),  A0(grid),
                        AB(grid),  A1(grid),
                                   A2(grid);

    DistMatrix<T,MC,MR> CT(grid),  C0(grid),
                        CB(grid),  C1(grid),
                                   C2(grid);

    // Temporary distributions
    DistMatrix<T,Star,MR> A1_Star_MR(grid);
    DistMatrix<T,Star,MC> D1_Star_MC(grid);
    DistMatrix<T,MR,  MC> D1_MR_MC(grid);
    DistMatrix<T,MC,  MR> D1(grid);

    // Start the algorithm
    BLAS::Scal( beta, C );
    LockedPartitionDown( A, AT,
                            AB );
    PartitionDown( C, CT,
                      CB );
    while( AB.Height() > 0 )
    {
        LockedRepartitionDown( AT,  A0,
                              /**/ /**/
                                    A1,
                               AB,  A2 );

        RepartitionDown( CT,  C0,
                        /**/ /**/
                              C1,
                         CB,  C2 );

        A1_Star_MR.ConformWith( B );
        D1_Star_MC.AlignWith( B );
        D1_Star_MC.ResizeTo( C1.Height(), C1.Width() );
        D1.AlignWith( C1 );
        //--------------------------------------------------------------------//
        A1_Star_MR = A1; // A1[*,MR] <- A1[MC,MR]

        // D1[*,MC] := alpha A1[*,MR] (B[MC,MR])^T
        //           = alpha A1[*,MR] (B^T)[MR,MC]
        BLAS::Gemm( Normal, orientationOfB, 
                    alpha, A1_Star_MR.LockedLocalMatrix(),
                           B.LockedLocalMatrix(),
                    (T)0,  D1_Star_MC.LocalMatrix()       );

        // C1[MC,MR] += scattered & transposed D1[*,MC] summed over grid rows
        D1_MR_MC.ReduceScatterFrom( D1_Star_MC );
        D1 = D1_MR_MC; 
        BLAS::Axpy( (T)1, D1, C1 );
        //--------------------------------------------------------------------//
        A1_Star_MR.FreeConstraints();
        D1_Star_MC.FreeConstraints();
        D1.FreeConstraints();

        SlideLockedPartitionDown( AT,  A0,
                                       A1,
                                 /**/ /**/
                                  AB,  A2 );
    
        SlidePartitionDown( CT,  C0,
                                 C1,
                           /**/ /**/
                            CB,  C2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Normal Transpose Gemm that avoids communicating the matrix C.
template<typename T>
void
Elemental::BLAS::Internal::GemmNTC
( const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::GemmNTC");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "{A,B,C} must be distributed over the same grid.";
    if( orientationOfB == Normal )
        throw "GemmNTC requires that B be (Conjugate)Transposed.";
    if( A.Height() != C.Height() ||
        B.Height() != C.Width()  ||
        A.Width()  != B.Width()    )
    {
        ostringstream msg;
        msg << "Nonconformal GemmNTC: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B ~ " << B.Height() << " x " << B.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        throw msg.str();
    }
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> AL(grid), AR(grid),
                        A0(grid), A1(grid), A2(grid);

    DistMatrix<T,MC,MR> BL(grid), BR(grid),
                        B0(grid), B1(grid), B2(grid);

    // Temporary distributions
    DistMatrix<T,MC,Star> A1_MC_Star(grid);
    DistMatrix<T,MR,Star> B1_MR_Star(grid);

    // Start the algorithm
    BLAS::Scal( beta, C );
    LockedPartitionRight( A, AL, AR );
    LockedPartitionRight( B, BL, BR );
    while( AR.Width() > 0 )
    {
        LockedRepartitionRight( AL, /**/ AR,
                                A0, /**/ A1, A2 );    

        LockedRepartitionRight( BL, /**/ BR,
                                B0, /**/ B1, B2 );

        A1_MC_Star.AlignWith( C );
        B1_MR_Star.AlignWith( C );
        //--------------------------------------------------------------------//
        A1_MC_Star = A1; // A1[MC,*] <- A1[MC,MR]
        B1_MR_Star = B1; // B1[MR,*] <- B1[MC,MR]

        // C[MC,MR] += alpha A1[MC,*] (B1[MR,*])^T
        //           = alpha A1[MC,*] (B1^T)[*,MR]
        BLAS::Gemm( Normal, orientationOfB, 
                    alpha, A1_MC_Star.LockedLocalMatrix(),
                           B1_MR_Star.LockedLocalMatrix(),
                    (T)1,  C.LocalMatrix()                );
        //--------------------------------------------------------------------//
        A1_MC_Star.FreeConstraints();
        B1_MR_Star.FreeConstraints();
 
        SlideLockedPartitionRight( AL,     /**/ AR,
                                   A0, A1, /**/ A2 );

        SlideLockedPartitionRight( BL,     /**/ BR,
                                   B0, B1, /**/ B2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Internal::GemmNT
( const Orientation orientationOfB,
  const float alpha, const DistMatrix<float,MC,MR>& A,            
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmNT
( const Orientation orientationOfB,
  const double alpha, const DistMatrix<double,MC,MR>& A,         
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Internal::GemmNT
( const Orientation orientationOfB,
  const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,  
                        const DistMatrix<scomplex,MC,MR>& B,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmNT
( const Orientation orientationOfB,
  const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,  
                        const DistMatrix<dcomplex,MC,MR>& B,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif
