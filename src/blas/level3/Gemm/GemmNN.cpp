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
elemental::blas::internal::GemmNN
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmNN");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "{A,B,C} must be distributed over the same grid.";
#endif
    const int m = C.Height();
    const int n = C.Width();
    const int k = A.Width();
    const float weightTowardsC = 2.0;
    const float weightAwayFromDot = 10.0;

    if( weightAwayFromDot*m <= k && weightAwayFromDot*n <= k )
    {
#ifndef RELEASE
        if( A.GetGrid().VCRank() == 0 )
            cout << "  GemmNN routing to GemmNNDot." << endl;
#endif
        blas::internal::GemmNNDot( alpha, A, B, beta, C );
    }
    else if( m <= n && weightTowardsC*m <= k )
    {
#ifndef RELEASE
        if( A.GetGrid().VCRank() == 0 )
            cout << "  GemmNN routing to GemmNNB." << endl;
#endif
        blas::internal::GemmNNB( alpha, A, B, beta, C );    
    }
    else if( n <= m && weightTowardsC*n <= k )
    {
#ifndef RELEASE
        if( A.GetGrid().VCRank() == 0 )
            cout << "  GemmNN routing to GemmNNA." << endl;
#endif
        blas::internal::GemmNNA( alpha, A, B, beta, C );
    }
    else
    {
#ifndef RELEASE
        if( A.GetGrid().VCRank() == 0 )
            cout << "  GemmNN routing to GemmNNC." << endl;
#endif
        blas::internal::GemmNNC( alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Normal Normal Gemm that avoids communicating the matrix A.
template<typename T>
void
elemental::blas::internal::GemmNNA
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmNNA");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "{A,B,C} must be distributed over the same grid.";
    if( A.Height() != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Width()  != B.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal GemmNNA: " << endl
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
    DistMatrix<T,VR,Star> B1_VR_Star(grid);
    DistMatrix<T,Star,MR> B1Trans_Star_MR(grid);
    DistMatrix<T,MC,Star> D1_MC_Star(grid);

    // Start the algorithm
    blas::Scal( beta, C );
    LockedPartitionRight( B, BL, BR );
    PartitionRight( C, CL, CR );
    while( BR.Width() > 0 )
    {
        LockedRepartitionRight
        ( BL, /**/     BR,
          B0, /**/ B1, B2 );

        RepartitionRight
        ( CL, /**/     CR,
          C0, /**/ C1, C2 );

        B1_VR_Star.AlignWith( A );
        B1Trans_Star_MR.AlignWith( A );
        D1_MC_Star.AlignWith( A );
        D1_MC_Star.ResizeTo( C1.Height(), C1.Width() );
        //--------------------------------------------------------------------//
        B1_VR_Star = B1;
        B1Trans_Star_MR.TransposeFrom( B1_VR_Star );

        // D1[MC,*] := alpha A[MC,MR] B1[MR,*]
        blas::Gemm
        ( Normal, Transpose, 
          alpha, A.LockedLocalMatrix(),
                 B1Trans_Star_MR.LockedLocalMatrix(),
          (T)0,  D1_MC_Star.LocalMatrix() );

        // C1[MC,MR] += scattered result of D1[MC,*] summed over grid rows
        C1.SumScatterUpdate( (T)1, D1_MC_Star );
        //--------------------------------------------------------------------//
        B1_VR_Star.FreeAlignments();
        B1Trans_Star_MR.FreeAlignments();
        D1_MC_Star.FreeAlignments();

        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );

        SlidePartitionRight
        ( CL,     /**/ CR,
          C0, C1, /**/ C2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Normal Normal Gemm that avoids communicating the matrix B.
template<typename T>
void 
elemental::blas::internal::GemmNNB
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmNNB");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "{A,B,C} must be distributed over the same grid.";
    if( A.Height() != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Width()  != B.Height()   )
    {
        ostringstream msg;
        msg << "Nonconformal GemmNNB: " << endl
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
    DistMatrix<T,Star,MC> A1_Star_MC(grid);
    DistMatrix<T,Star,MR> D1_Star_MR(grid);

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

        A1_Star_MC.AlignWith( B );
        D1_Star_MR.AlignWith( B );
        D1_Star_MR.ResizeTo( C1.Height(), C1.Width() );
        //--------------------------------------------------------------------//
        A1_Star_MC = A1; // A1[*,MC] <- A1[MC,MR]

        // D1[*,MR] := alpha A1[*,MC] B[MC,MR]
        blas::Gemm
        ( Normal, Normal, 
          alpha, A1_Star_MC.LockedLocalMatrix(),
                 B.LockedLocalMatrix(),
          (T)0,  D1_Star_MR.LocalMatrix() );

        // C1[MC,MR] += scattered result of D1[*,MR] summed over grid cols
        C1.SumScatterUpdate( (T)1, D1_Star_MR );
        //--------------------------------------------------------------------//
        A1_Star_MC.FreeAlignments();
        D1_Star_MR.FreeAlignments();

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

// Normal Normal Gemm that avoids communicating the matrix C.
template<typename T>
void 
elemental::blas::internal::GemmNNC
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmNNC");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "{A,B,C} must be distributed over the same grid.";
    if( A.Height() != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Width()  != B.Height()   )
    {
        ostringstream msg;
        msg << "Nonconformal GemmNNC: " << endl
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

    DistMatrix<T,MC,MR> BT(grid),  B0(grid),
                        BB(grid),  B1(grid),
                                   B2(grid);

    // Temporary distributions
    DistMatrix<T,MC,Star> A1_MC_Star(grid);
    DistMatrix<T,MR,Star> B1Trans_MR_Star(grid); 

    // Start the algorithm
    blas::Scal( beta, C );
    LockedPartitionRight( A, AL, AR ); 
    LockedPartitionDown
    ( B, BT, 
         BB ); 
    while( AR.Width() > 0 )
    {
        LockedRepartitionRight( AL, /**/ AR,
                                A0, /**/ A1, A2 );

        LockedRepartitionDown( BT,  B0,
                              /**/ /**/
                                    B1, 
                               BB,  B2 );

        A1_MC_Star.AlignWith( C );
        B1Trans_MR_Star.AlignWith( C );
        //--------------------------------------------------------------------//
        A1_MC_Star = A1; 
        B1Trans_MR_Star.TransposeFrom( B1 );

        // C[MC,MR] += alpha A1[MC,*] (B1^T[MR,*])^T
        //           = alpha A1[MC,*] B1[*,MR]
        blas::Gemm
        ( Normal, Transpose, 
          alpha, A1_MC_Star.LockedLocalMatrix(),
                 B1Trans_MR_Star.LockedLocalMatrix(),
          (T)1,  C.LocalMatrix() );
        //--------------------------------------------------------------------//
        A1_MC_Star.FreeAlignments();
        B1Trans_MR_Star.FreeAlignments();

        SlideLockedPartitionRight( AL,     /**/ AR,
                                   A0, A1, /**/ A2 );

        SlideLockedPartitionDown( BT,  B0,
                                       B1,
                                 /**/ /**/
                                  BB,  B2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Normal Normal Gemm for panel-panel dot products. 
template<typename T>
void 
elemental::blas::internal::GemmNNDot
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmNNDot");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw "{A,B,C} must be distributed over the same grid.";
    if( A.Height() != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Width()  != B.Height()   )
    {
        ostringstream msg;
        msg << "Nonconformal GemmNNDot: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B ~ " << B.Height() << " x " << B.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width();
        const string& s = msg.str();
        throw s.c_str();
    }
#endif
    const Grid& grid = A.GetGrid();

    if( A.Height() > B.Width() )
    {
        // Matrix views
        DistMatrix<T,MC,MR> AT(grid), AB(grid),
                            A0(grid), A1(grid), A2(grid);         

        DistMatrix<T,MC,MR> BL(grid),  B0(grid),
                            BR(grid),  B1(grid),
                                       B2(grid);

        DistMatrix<T,MC,MR> CT(grid), C0(grid), C1L(grid), C1R(grid),
                            CB(grid), C1(grid), C10(grid), C11(grid), C12(grid),
                                      C2(grid);

        // Temporary distributions
        DistMatrix<T,Star,VC> A1_Star_VC(grid);
        DistMatrix<T,VC,Star> B1_VC_Star(grid);
        DistMatrix<T,Star,Star> C11_Star_Star(grid);

        // Star the algorithm
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

            A1_Star_VC.AlignWith( B1 );
            //----------------------------------------------------------------//
            A1_Star_VC = A1; 
            //----------------------------------------------------------------//

            LockedPartitionRight( B, BL, BR );
            PartitionRight( C1, C1L, C1R );
            while( BR.Width() > 0 )
            {
                LockedRepartitionRight
                ( BL, /**/ BR,
                  B0, /**/ B1, B2 );

                RepartitionRight
                ( C1L, /**/ C1R,
                  C10, /**/ C11, C12 );

                B1_VC_Star.AlignWith( A1_Star_VC );
                C11_Star_Star.ResizeTo( C11.Height(), C11.Width() );
                //------------------------------------------------------------//
                B1_VC_Star = B1;
                blas::Gemm
                ( Normal, Normal,
                  alpha, A1_Star_VC.LockedLocalMatrix(),
                         B1_VC_Star.LockedLocalMatrix(),
                  (T)0,  C11_Star_Star.LocalMatrix() );
                C11.SumScatterUpdate( (T)1, C11_Star_Star );
                //------------------------------------------------------------//
                B1_VC_Star.FreeAlignments();

                SlideLockedPartitionRight
                ( BL,     /**/ BR,
                  B0, B1, /**/ B2 );

                SlidePartitionRight
                ( C1L,      /**/ C1R,
                  C10, C11, /**/ C12 );
            }
            A1_Star_VC.FreeAlignments();

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
    }
    else
    {
        // Matrix views
        DistMatrix<T,MC,MR> AT(grid), AB(grid),
                            A0(grid), A1(grid), A2(grid);         

        DistMatrix<T,MC,MR> BL(grid),  B0(grid),
                            BR(grid),  B1(grid),
                                       B2(grid);

        DistMatrix<T,MC,MR> 
            CL(grid), CR(grid),            C1T(grid),  C01(grid),
            C0(grid), C1(grid), C2(grid),  C1B(grid),  C11(grid),
                                                       C21(grid);

        // Temporary distributions
        DistMatrix<T,Star,VR> A1_Star_VR(grid);
        DistMatrix<T,VR,Star> B1_VR_Star(grid);
        DistMatrix<T,Star,Star> C11_Star_Star(grid);

        // Star the algorithm
        blas::Scal( beta, C );
        LockedPartitionRight( B, BL, BR );
        PartitionRight( C, CL, CR );
        while( BR.Width() > 0 )
        {
            LockedRepartitionRight
            ( BL, /**/ BR,
              B0, /**/ B1, B2 );

            RepartitionRight
            ( CL, /**/ CR,
              C0, /**/ C1, C2 );

            B1_VR_Star.AlignWith( A1 );
            //----------------------------------------------------------------//
            B1_VR_Star = B1;
            //----------------------------------------------------------------//

            LockedPartitionDown
            ( A, AT,
                 AB );
            PartitionDown
            ( C1, C1T,
                  C1B );
            while( AB.Height() > 0 )
            {
                LockedRepartitionDown
                ( AT,  A0,
                 /**/ /**/
                       A1,
                  AB,  A2 );

                RepartitionDown
                ( C1T,  C01,
                 /***/ /***/
                        C11,
                  C1B,  C21 );

                A1_Star_VR.AlignWith( B1_VR_Star );
                C11_Star_Star.ResizeTo( C11.Height(), C11.Width() );
                //------------------------------------------------------------//
                A1_Star_VR = A1;
                blas::Gemm
                ( Normal, Normal,
                  alpha, A1_Star_VR.LockedLocalMatrix(),
                         B1_VR_Star.LockedLocalMatrix(),
                  (T)0,  C11_Star_Star.LocalMatrix() );
                C11.SumScatterUpdate( (T)1, C11_Star_Star );
                //------------------------------------------------------------//
                A1_Star_VR.FreeAlignments();

                SlideLockedPartitionDown
                ( AT,  A0,
                       A1,
                 /**/ /**/
                  AB,  A2 );

                SlidePartitionDown
                ( C1T,  C01,
                        C11,
                 /***/ /***/
                  C1B,  C21 );
            }
            B1_VR_Star.FreeAlignments();

            SlideLockedPartitionRight
            ( BL,     /**/ BR,
              B0, B1, /**/ B2 ); 

            SlidePartitionRight
            ( CL,     /**/ CR,
              C0, C1, /**/ C2 );
        }
    }

#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::GemmNN
( float alpha, const DistMatrix<float,MC,MR>& A,     
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::GemmNN
( double alpha, const DistMatrix<double,MC,MR>& A,         
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::GemmNN
( scomplex alpha, const DistMatrix<scomplex,MC,MR>& A, 
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::GemmNN
( dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A, 
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

