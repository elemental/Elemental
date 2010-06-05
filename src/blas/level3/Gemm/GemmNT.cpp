/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::blas::internal::GemmNT
( Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmNT");
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
        blas::internal::GemmNTB( orientationOfB, alpha, A, B, beta, C );
    }
    else if( n <= m && weightTowardsC*n <= k )
    {
#ifndef RELEASE
        if( A.GetGrid().VCRank() == 0 )
            cout << "  GemmNT routing to GemmNTA." << endl;
#endif
        blas::internal::GemmNTA( orientationOfB, alpha, A, B, beta, C );
    }
    else
    {
#ifndef RELEASE
        if( A.GetGrid().VCRank() == 0 )
            cout << "  GemmNT routing to GemmNTC." << endl;
#endif
        blas::internal::GemmNTC( orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Normal Transpose Gemm that avoids communicating the matrix A.
template<typename T>
void
elemental::blas::internal::GemmNTA
( Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmNTA");
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
        const string& s = msg.str();
        throw s.c_str();
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
    blas::Scal( beta, C );
    LockedPartitionDown
    ( B, BT,
         BB );
    PartitionRight( C, CL, CR );
    while( BB.Height() > 0 )
    {
        LockedRepartitionDown
        ( BT,  B0,
         /**/ /**/
               B1,
          BB,  B2 );

        RepartitionRight
        ( CL, /**/     CR,
          C0, /**/ C1, C2 );

        B1_Star_MR.AlignWith( A );
        D1_MC_Star.AlignWith( A );
        D1_MC_Star.ResizeTo( C1.Height(), C1.Width() );
        //--------------------------------------------------------------------//
        B1_Star_MR = B1; // B1[*,MR] <- B1[MC,MR]

        // C1[MC,*] := alpha A[MC,MR] (B1[*,MR])^T
        //           = alpha A[MC,MR] (B1^T)[MR,*]
        blas::Gemm
        ( Normal, orientationOfB, 
          alpha, A.LockedLocalMatrix(),
                 B1_Star_MR.LockedLocalMatrix(),
          (T)0,  D1_MC_Star.LocalMatrix() );

        // C1[MC,MR] += scattered result of D1[MC,*] summed over grid rows
        C1.SumScatterUpdate( (T)1, D1_MC_Star );
        //--------------------------------------------------------------------//
        B1_Star_MR.FreeAlignments();
        D1_MC_Star.FreeAlignments();

        SlideLockedPartitionDown
        ( BT,  B0,
               B1,
         /**/ /**/
          BB,  B2 );

        SlidePartitionRight
        ( CL,     /**/ CR,
          C0, C1, /**/ C2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Normal Transpose Gemm that avoids communicating the matrix B.
template<typename T>
void
elemental::blas::internal::GemmNTB
( Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmNTB");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "{A,B,C} must be distributed over the same grid.";
    if( orientationOfB == Normal )
        throw "GemmNTB requires that B be (Conjugate)Transposed.";
    if( A.Height() != C.Height() ||
        B.Height() != C.Width()  ||
        A.Width()  != B.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal GemmNTB: " << endl
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

    DistMatrix<T,MC,MR> CT(grid),  C0(grid),
                        CB(grid),  C1(grid),
                                   C2(grid);

    // Temporary distributions
    DistMatrix<T,Star,MR> A1_Star_MR(grid);
    DistMatrix<T,Star,MC> D1_Star_MC(grid);
    DistMatrix<T,MR,  MC> D1_MR_MC(grid);
    DistMatrix<T,MC,  MR> D1(grid);

    // Start the algorithm
    blas::Scal( beta, C );
    LockedPartitionDown
    ( A, AT,
         AB );
    PartitionDown
    ( C, CT,
         CB );
    while( AB.Height() > 0 )
    {
        LockedRepartitionDown
        ( AT,  A0,
         /**/ /**/
               A1,
          AB,  A2 );

        RepartitionDown
        ( CT,  C0,
         /**/ /**/
               C1,
          CB,  C2 );

        A1_Star_MR.AlignWith( B );
        D1_Star_MC.AlignWith( B );
        D1_Star_MC.ResizeTo( C1.Height(), C1.Width() );
        D1.AlignWith( C1 );
        //--------------------------------------------------------------------//
        A1_Star_MR = A1; // A1[*,MR] <- A1[MC,MR]

        // D1[*,MC] := alpha A1[*,MR] (B[MC,MR])^T
        //           = alpha A1[*,MR] (B^T)[MR,MC]
        blas::Gemm
        ( Normal, orientationOfB, 
          alpha, A1_Star_MR.LockedLocalMatrix(),
                 B.LockedLocalMatrix(),
          (T)0,  D1_Star_MC.LocalMatrix() );

        // C1[MC,MR] += scattered & transposed D1[*,MC] summed over grid rows
        D1_MR_MC.SumScatterFrom( D1_Star_MC );
        D1 = D1_MR_MC; 
        blas::Axpy( (T)1, D1, C1 );
        //--------------------------------------------------------------------//
        A1_Star_MR.FreeAlignments();
        D1_Star_MC.FreeAlignments();
        D1.FreeAlignments();

        SlideLockedPartitionDown
        ( AT,  A0,
               A1,
         /**/ /**/
          AB,  A2 );
    
        SlidePartitionDown
        ( CT,  C0,
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
elemental::blas::internal::GemmNTC
( Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmNTC");
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
        const string& s = msg.str();
        throw s.c_str();
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
    blas::Scal( beta, C );
    LockedPartitionRight( A, AL, AR );
    LockedPartitionRight( B, BL, BR );
    while( AR.Width() > 0 )
    {
        LockedRepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );    

        LockedRepartitionRight
        ( BL, /**/ BR,
          B0, /**/ B1, B2 );

        A1_MC_Star.AlignWith( C );
        B1_MR_Star.AlignWith( C );
        //--------------------------------------------------------------------//
        A1_MC_Star = A1; // A1[MC,*] <- A1[MC,MR]
        B1_MR_Star = B1; // B1[MR,*] <- B1[MC,MR]

        // C[MC,MR] += alpha A1[MC,*] (B1[MR,*])^T
        //           = alpha A1[MC,*] (B1^T)[*,MR]
        blas::Gemm
        ( Normal, orientationOfB, 
          alpha, A1_MC_Star.LockedLocalMatrix(),
                 B1_MR_Star.LockedLocalMatrix(),
          (T)1,  C.LocalMatrix() );
        //--------------------------------------------------------------------//
        A1_MC_Star.FreeAlignments();
        B1_MR_Star.FreeAlignments();
 
        SlideLockedPartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );

        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::GemmNT
( Orientation orientationOfB,
  float alpha, const DistMatrix<float,MC,MR>& A,            
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::GemmNT
( Orientation orientationOfB,
  double alpha, const DistMatrix<double,MC,MR>& A,         
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::GemmNT
( Orientation orientationOfB,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,  
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::GemmNT
( Orientation orientationOfB,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,  
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

