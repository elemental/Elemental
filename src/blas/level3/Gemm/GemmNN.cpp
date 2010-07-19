/*
   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Elemental.  If not, see <http://www.gnu.org/licenses/>.
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
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
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
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
    if( A.Height() != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Width()  != B.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal GemmNNA: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B ~ " << B.Height() << " x " << B.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> BL(g), BR(g),
                        B0(g), B1(g), B2(g);
    DistMatrix<T,MC,MR> CL(g), CR(g),
                        C0(g), C1(g), C2(g);

    // Temporary distributions
    DistMatrix<T,VR,Star> B1_VR_Star(g);
    DistMatrix<T,Star,MR> B1Trans_Star_MR(g);
    DistMatrix<T,MC,Star> D1_MC_Star(g);

    // Start the algorithm
    blas::Scal( beta, C );
    LockedPartitionRight( B, BL, BR, 0 );
    PartitionRight( C, CL, CR, 0 );
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
        blas::internal::LocalGemm
        ( Normal, Transpose, alpha, A, B1Trans_Star_MR, (T)0, D1_MC_Star );

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
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
    if( A.Height() != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Width()  != B.Height()   )
    {
        ostringstream msg;
        msg << "Nonconformal GemmNNB: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B ~ " << B.Height() << " x " << B.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> AT(g),  A0(g),
                        AB(g),  A1(g),
                                A2(g);
    DistMatrix<T,MC,MR> CT(g),  C0(g),
                        CB(g),  C1(g),
                                C2(g);

    // Temporary distributions
    DistMatrix<T,Star,MC> A1_Star_MC(g);
    DistMatrix<T,Star,MR> D1_Star_MR(g);

    // Start the algorithm
    blas::Scal( beta, C );
    LockedPartitionDown
    ( A, AT,
         AB, 0 );
    PartitionDown
    ( C, CT,
         CB, 0 );
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
        blas::internal::LocalGemm
        ( Normal, Normal, alpha, A1_Star_MC, B, (T)0, D1_Star_MR );

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
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
    if( A.Height() != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Width()  != B.Height()   )
    {
        ostringstream msg;
        msg << "Nonconformal GemmNNC: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B ~ " << B.Height() << " x " << B.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> AL(g), AR(g),
                        A0(g), A1(g), A2(g);         

    DistMatrix<T,MC,MR> BT(g),  B0(g),
                        BB(g),  B1(g),
                                B2(g);

    // Temporary distributions
    DistMatrix<T,MC,Star> A1_MC_Star(g);
    DistMatrix<T,MR,Star> B1Trans_MR_Star(g); 

    // Start the algorithm
    blas::Scal( beta, C );
    LockedPartitionRight( A, AL, AR, 0 ); 
    LockedPartitionDown
    ( B, BT, 
         BB, 0 ); 
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
        blas::internal::LocalGemm
        ( Normal, Transpose, alpha, A1_MC_Star, B1Trans_MR_Star, (T)1, C );
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
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
    if( A.Height() != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Width()  != B.Height()   )
    {
        ostringstream msg;
        msg << "Nonconformal GemmNNDot: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B ~ " << B.Height() << " x " << B.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width();
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = A.GetGrid();

    if( A.Height() > B.Width() )
    {
        // Matrix views
        DistMatrix<T,MC,MR> AT(g), AB(g),
                            A0(g), A1(g), A2(g);         

        DistMatrix<T,MC,MR> BL(g),  B0(g),
                            BR(g),  B1(g),
                                    B2(g);

        DistMatrix<T,MC,MR> CT(g), C0(g), C1L(g), C1R(g),
                            CB(g), C1(g), C10(g), C11(g), C12(g),
                                   C2(g);

        // Temporary distributions
        DistMatrix<T,Star,VC> A1_Star_VC(g);
        DistMatrix<T,VC,Star> B1_VC_Star(g);
        DistMatrix<T,Star,Star> C11_Star_Star(g);

        // Star the algorithm
        blas::Scal( beta, C );
        LockedPartitionDown
        ( A, AT,
             AB, 0 );
        PartitionDown
        ( C, CT,
             CB, 0 );
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

            LockedPartitionRight( B, BL, BR, 0 );
            PartitionRight( C1, C1L, C1R, 0 );
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
                blas::internal::LocalGemm
                ( Normal, Normal, 
                  alpha, A1_Star_VC, B1_VC_Star, (T)0, C11_Star_Star );
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
        DistMatrix<T,MC,MR> AT(g), AB(g),
                            A0(g), A1(g), A2(g);         

        DistMatrix<T,MC,MR> BL(g),  B0(g),
                            BR(g),  B1(g),
                                    B2(g);

        DistMatrix<T,MC,MR> 
            CL(g), CR(g),         C1T(g),  C01(g),
            C0(g), C1(g), C2(g),  C1B(g),  C11(g),
                                           C21(g);

        // Temporary distributions
        DistMatrix<T,Star,VR> A1_Star_VR(g);
        DistMatrix<T,VR,Star> B1_VR_Star(g);
        DistMatrix<T,Star,Star> C11_Star_Star(g);

        // Star the algorithm
        blas::Scal( beta, C );
        LockedPartitionRight( B, BL, BR, 0 );
        PartitionRight( C, CL, CR, 0 );
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
                 AB, 0 );
            PartitionDown
            ( C1, C1T,
                  C1B, 0 );
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
                blas::internal::LocalGemm
                ( Normal, Normal,
                  alpha, A1_Star_VR, B1_VR_Star, (T)0, C11_Star_Star );
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

