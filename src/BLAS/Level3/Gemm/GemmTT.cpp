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
Elemental::BLAS::Internal::GemmTT
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::GemmTT");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
    {
        if( A.GetGrid().VCRank() == 0 )
            cerr << "{A,B,C} must be distributed over the same grid." << endl;
        DumpCallStack();
        throw exception();
    }
    if( orientationOfA == Normal || orientationOfB == Normal )
    {
        if( A.GetGrid().VCRank() == 0 )
            cerr << "GemmTT expects A and B to be transposed." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    const int m = C.Height();
    const int n = C.Width();
    const int k = A.Height();
    const float weightTowardsC = 2.0;

    if( m <= n && weightTowardsC*m <= k )
    {
#ifndef RELEASE
        if( A.GetGrid().VCRank() == 0 )
            cout << "  GemmTT routing to GemmTTB." << endl;
#endif
        BLAS::Internal::GemmTTB( orientationOfA, orientationOfB, 
                                  alpha, A, B, beta, C );
    }
    else if( n <= m && weightTowardsC*n <= k )
    {
#ifndef RELEASE
        if( A.GetGrid().VCRank() == 0 )
            cout << "  GemmTT routing to GemmTTA." << endl;
#endif
        BLAS::Internal::GemmTTA( orientationOfA, orientationOfB,
                                  alpha, A, B, beta, C );
    }
    else
    {
#ifndef RELEASE
        if( A.GetGrid().VCRank() == 0 )
            cout << "  GemmTT routing to GemmTTC." << endl;
#endif
        BLAS::Internal::GemmTTC( orientationOfA, orientationOfB,
                                  alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Transpose Transpose Gemm that avoids communicating the matrix A.
template<typename T>
void
Elemental::BLAS::Internal::GemmTTA
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::GemmTTA");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
    {
        if( A.GetGrid().VCRank() == 0 )
            cerr << "{A,B,C} must be distributed over the same grid." << endl;
        DumpCallStack();
        throw exception();
    }
    if( orientationOfA == Normal || orientationOfB == Normal )
    {
        if( A.GetGrid().VCRank() == 0 )
        {
            cerr << "GemmTTA expects A and B to be (Conjugate)Transposed." 
                 << endl;
        }
        DumpCallStack();
        throw exception();
    }
    if( A.Width()  != C.Height() ||
        B.Height() != C.Width()  ||
        A.Height() != B.Width()    )
    {
        if( A.GetGrid().VCRank() == 0 )
        {
            cerr << "Nonconformal GemmTTA: " <<
            endl << "  A ~ " << A.Height() << " x " << A.Width() <<
            endl << "  B ~ " << B.Height() << " x " << B.Width() <<
            endl << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        }
        DumpCallStack();
        throw exception();
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
    DistMatrix<T,Star,MC  > B1_Star_MC(grid);
    DistMatrix<T,MR,  Star> D1_MR_Star(grid);
    DistMatrix<T,MR,  MC  > D1_MR_MC(grid);
    DistMatrix<T,MC,  MR  > D1(grid);

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

        B1_Star_MC.ConformWith( A ); 
        D1_MR_Star.AlignWith( A );  
        D1_MR_Star.ResizeTo( C1.Height(), C1.Width() );
        D1.AlignWith( C1 );  
        //--------------------------------------------------------------------//
        B1_Star_MC = B1; // B1[*,MC] <- B1[MC,MR]

        // D1[MR,*] := alpha (A[MC,MR])^T (B1[*,MC])^T
        //           = alpha (A^T)[MR,MC] (B1^T)[MC,*]
        BLAS::Gemm( orientationOfA, orientationOfB,
                    alpha, A.LockedLocalMatrix(),
                           B1_Star_MC.LockedLocalMatrix(),
                    (T)0,  D1_MR_Star.LocalMatrix()       );

        // C1[MC,MR] += scattered & transposed D1[MR,*] summed over grid cols
        D1_MR_MC.ReduceScatterFrom ( D1_MR_Star );
        D1 = D1_MR_MC; 
        BLAS::Axpy( (T)1, D1, C1 );
        //--------------------------------------------------------------------//
        B1_Star_MC.FreeConstraints();
        D1_MR_Star.FreeConstraints();
        D1.FreeConstraints();
        
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

// Transpose Transpose Gemm that avoids communicating the matrix B.
template<typename T>
void
Elemental::BLAS::Internal::GemmTTB
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::GemmTTB");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
    {
        if( A.GetGrid().VCRank() == 0 )
            cerr << "{A,B,C} must be distributed over the same grid." << endl;
        DumpCallStack();
        throw exception();
    }
    if( orientationOfA == Normal || orientationOfB == Normal )
    {
        if( A.GetGrid().VCRank() == 0 )
        {
            cerr << "GemmTTB expects A and B to be (Conjugate)Transposed." 
                 << endl;
        }
        DumpCallStack();
        throw exception();
    }
    if( A.Width()  != C.Height() ||
        B.Height() != C.Width()  ||
        A.Height() != B.Width()    )
    {
        if( A.GetGrid().VCRank() == 0 )
        {
            cerr << "Nonconformal GemmTTB: " <<
            endl << "  A ~ " << A.Height() << " x " << A.Width() <<
            endl << "  B ~ " << B.Height() << " x " << B.Width() <<
            endl << "  C ~ " << C.Height() << " x " << C.Width() <<
            endl;
        }
        DumpCallStack();
        throw exception();
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
    DistMatrix<T,MR,  Star> A1_MR_Star(grid);
    DistMatrix<T,Star,MC  > D1_Star_MC(grid);
    DistMatrix<T,MR,  MC  > D1_MR_MC(grid);
    DistMatrix<T,MC,  MR  > D1(grid);

    // Start the algorithm 
    BLAS::Scal( beta, C );
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

        A1_MR_Star.ConformWith( B );
        D1_Star_MC.AlignWith( B );
        D1_Star_MC.ResizeTo( C1.Height(), C1.Width() );
        D1.AlignWith( C1 );
        //--------------------------------------------------------------------//
        A1_MR_Star = A1; // A1[MR,*] <- A1[MC,MR]
 
        // D1[*,MC] := alpha (A1[MR,*])^T (B[MC,MR])^T
        //           = alpha (A1^T)[*,MR] (B^T)[MR,MC]
        BLAS::Gemm( orientationOfA, orientationOfB,
                    alpha, A1_MR_Star.LockedLocalMatrix(),
                           B.LockedLocalMatrix(),
                    (T)0,  D1_Star_MC.LocalMatrix()       );

        // C1[MC,MR] += scattered & transposed D1[*,MC] summed over grid rows
        D1_MR_MC.ReduceScatterFrom( D1_Star_MC );
        D1 = D1_MR_MC; 
        BLAS::Axpy( (T)1, D1, C1 );
        //--------------------------------------------------------------------//
        A1_MR_Star.FreeConstraints();
        D1_Star_MC.FreeConstraints();
        D1.FreeConstraints();

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

// Transpose Transpose Gemm that avoids communicating the matrix C.
template<typename T>
void
Elemental::BLAS::Internal::GemmTTC
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::GemmTTC");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
    {
        if( A.GetGrid().VCRank() == 0 )
            cerr << "{A,B,C} must be distributed over the same grid." << endl;
        DumpCallStack();
        throw exception();
    }
    if( orientationOfA == Normal || orientationOfB == Normal )
    {
        if( A.GetGrid().VCRank() == 0 )
        {
            cerr << "GemmTTC expects A and B to be (Conjugate)Transposed." 
                 << endl;
        }
        DumpCallStack();
        throw exception();
    }
    if( A.Width()  != C.Height() ||
        B.Height() != C.Width()  ||
        A.Height() != B.Width()    )
    {
        if( A.GetGrid().VCRank() == 0 )
        {
            cerr << "Nonconformal GemmTTC: " <<
            endl << "  A ~ " << A.Height() << " x " << A.Width() <<
            endl << "  B ~ " << B.Height() << " x " << B.Width() <<
            endl << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        }
        DumpCallStack();
        throw exception();
    }
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> AT(grid),  A0(grid),
                        AB(grid),  A1(grid),
                                   A2(grid);

    DistMatrix<T,MC,MR> BL(grid), BR(grid),
                        B0(grid), B1(grid), B2(grid);

    // Temporary distributions
    DistMatrix<T,Star,MC> A1_Star_MC(grid);
    DistMatrix<T,MR,Star> B1_MR_Star(grid);
    
    // Start the algorithm    
    BLAS::Scal( beta, C );
    LockedPartitionDown( A, AT,
                            AB ); 
    LockedPartitionRight( B, BL, BR );
    while( AB.Height() > 0 )
    {
        LockedRepartitionDown( AT,  A0,
                              /**/ /**/
                                    A1,
                               AB,  A2 );

        LockedRepartitionRight( BL, /**/     BR,
                                B0, /**/ B1, B2 );

        A1_Star_MC.AlignWith( C );
        B1_MR_Star.AlignWith( C );
        //--------------------------------------------------------------------//
        A1_Star_MC = A1; // A1[*,MC] <- A1[MC,MR]
        B1_MR_Star = B1; // B1[MR,*] <- B1[MC,MR]

        // C[MC,MR] += alpha (A1[*,MC])^T (B1[MR,*])^T
        //           = alpha (A1^T)[MC,*] (B1^T)[*,MR]
        BLAS::Gemm( orientationOfA, orientationOfB, 
                    alpha, A1_Star_MC.LockedLocalMatrix(),
                           B1_MR_Star.LockedLocalMatrix(),
                    (T)1,  C.LocalMatrix()                );
        //--------------------------------------------------------------------//
        A1_Star_MC.FreeConstraints();
        B1_MR_Star.FreeConstraints();

        SlideLockedPartitionDown( AT,  A0,
                                       A1,
                                 /**/ /**/
                                  AB,  A2 );

        SlideLockedPartitionRight( BL,     /**/ BR,
                                   B0, B1, /**/ B2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Internal::GemmTT
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const float alpha, const DistMatrix<float,MC,MR>& A, 
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmTT
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const double alpha, const DistMatrix<double,MC,MR>& A, 
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Internal::GemmTT
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A, 
                        const DistMatrix<scomplex,MC,MR>& B,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::GemmTT
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A, 
                        const DistMatrix<dcomplex,MC,MR>& B,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif
