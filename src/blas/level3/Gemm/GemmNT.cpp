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
elemental::blas::internal::GemmNT
( Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmNT");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
    if( orientationOfB == Normal )
        throw logic_error( "GemmNT requires that B be (Conjugate)Transposed." );
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
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
    if( orientationOfB == Normal )
        throw logic_error
        ( "GemmTNA requires that B be (Conjugate)Transposed." );
    if( A.Height() != C.Height() ||
        B.Height() != C.Width()  ||
        A.Width()  != B.Width()    )
    {
        ostringstream msg;
        msg << "Nonconformal GemmNTA: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  B ~ " << B.Height() << " x " << B.Width() << endl
            << "  C ~ " << C.Height() << " x " << C.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> BT(g),  B0(g),
                        BB(g),  B1(g),
                                B2(g);

    DistMatrix<T,MC,MR> CL(g), CR(g),
                        C0(g), C1(g), C2(g);

    // Temporary distributions
    DistMatrix<T,Star,MR> B1_Star_MR(g);
    DistMatrix<T,MC,Star> D1_MC_Star(g);

    // Start the algorithm
    blas::Scal( beta, C );
    LockedPartitionDown
    ( B, BT,
         BB, 0 );
    PartitionRight( C, CL, CR, 0 );
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
        blas::internal::LocalGemm
        ( Normal, orientationOfB, alpha, A, B1_Star_MR, (T)0, D1_MC_Star );

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
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
    if( orientationOfB == Normal )
        throw logic_error
        ( "GemmNTB requires that B be (Conjugate)Transposed." );
    if( A.Height() != C.Height() ||
        B.Height() != C.Width()  ||
        A.Width()  != B.Width() )
    {
        ostringstream msg;
        msg << "Nonconformal GemmNTB: " << endl
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
    DistMatrix<T,Star,MR> A1_Star_MR(g);
    DistMatrix<T,Star,MC> D1_Star_MC(g);
    DistMatrix<T,MR,  MC> D1_MR_MC(g);
    DistMatrix<T,MC,  MR> D1(g);

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

        A1_Star_MR.AlignWith( B );
        D1_Star_MC.AlignWith( B );
        D1_Star_MC.ResizeTo( C1.Height(), C1.Width() );
        D1.AlignWith( C1 );
        //--------------------------------------------------------------------//
        A1_Star_MR = A1; // A1[*,MR] <- A1[MC,MR]

        // D1[*,MC] := alpha A1[*,MR] (B[MC,MR])^T
        //           = alpha A1[*,MR] (B^T)[MR,MC]
        blas::internal::LocalGemm
        ( Normal, orientationOfB, alpha, A1_Star_MR, B, (T)0, D1_Star_MC );

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
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
    if( orientationOfB == Normal )
        throw logic_error
        ( "GemmNTC requires that B be (Conjugate)Transposed." );
    if( A.Height() != C.Height() ||
        B.Height() != C.Width()  ||
        A.Width()  != B.Width()    )
    {
        ostringstream msg;
        msg << "Nonconformal GemmNTC: " << endl
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

    DistMatrix<T,MC,MR> BL(g), BR(g),
                        B0(g), B1(g), B2(g);

    // Temporary distributions
    DistMatrix<T,MC,Star> A1_MC_Star(g);
    DistMatrix<T,MR,Star> B1_MR_Star(g);

    // Start the algorithm
    blas::Scal( beta, C );
    LockedPartitionRight( A, AL, AR, 0 );
    LockedPartitionRight( B, BL, BR, 0 );
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
        blas::internal::LocalGemm
        ( Normal, orientationOfB, alpha, A1_MC_Star, B1_MR_Star, (T)1, C );
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

