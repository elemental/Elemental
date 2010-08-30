/*
   Copyright (c) 2009-2010, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::blas::internal::GemmTT
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmTT");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
    if( orientationOfA == Normal || orientationOfB == Normal )
        throw logic_error( "GemmTT expects A and B to be transposed." );
#endif
    const int m = C.Height();
    const int n = C.Width();
    const int k = A.Height();
    const float weightTowardsC = 2.0;

    if( m <= n && weightTowardsC*m <= k )
    {
        blas::internal::GemmTTB
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
    else if( n <= m && weightTowardsC*n <= k )
    {
        blas::internal::GemmTTA
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
    else
    {
        blas::internal::GemmTTC
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Transpose Transpose Gemm that avoids communicating the matrix A.
template<typename T>
void
elemental::blas::internal::GemmTTA
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmTTA");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
    if( orientationOfA == Normal || orientationOfB == Normal )
        throw logic_error
        ( "GemmTTA expects A and B to be (Conjugate)Transposed." );
    if( A.Width()  != C.Height() ||
        B.Height() != C.Width()  ||
        A.Height() != B.Width()    )
    {
        ostringstream msg;
        msg << "Nonconformal GemmTTA: " << endl
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
    DistMatrix<T,Star,MC  > B1_Star_MC(g);
    DistMatrix<T,MR,  Star> D1_MR_Star(g);
    DistMatrix<T,MR,  MC  > D1_MR_MC(g);
    DistMatrix<T,MC,  MR  > D1(g);

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

        B1_Star_MC.AlignWith( A ); 
        D1_MR_Star.AlignWith( A );  
        D1_MR_Star.ResizeTo( C1.Height(), C1.Width() );
        D1.AlignWith( C1 );  
        //--------------------------------------------------------------------//
        B1_Star_MC = B1; // B1[*,MC] <- B1[MC,MR]

        // D1[MR,*] := alpha (A[MC,MR])^T (B1[*,MC])^T
        //           = alpha (A^T)[MR,MC] (B1^T)[MC,*]
        blas::internal::LocalGemm
        ( orientationOfA, orientationOfB, 
          alpha, A, B1_Star_MC, (T)0, D1_MR_Star );

        // C1[MC,MR] += scattered & transposed D1[MR,*] summed over grid cols
        D1_MR_MC.SumScatterFrom( D1_MR_Star );
        D1 = D1_MR_MC; 
        blas::Axpy( (T)1, D1, C1 );
        //--------------------------------------------------------------------//
        B1_Star_MC.FreeAlignments();
        D1_MR_Star.FreeAlignments();
        D1.FreeAlignments();
        
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

// Transpose Transpose Gemm that avoids communicating the matrix B.
template<typename T>
void
elemental::blas::internal::GemmTTB
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmTTB");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
    if( orientationOfA == Normal || orientationOfB == Normal )
        throw logic_error
        ( "GemmTTB expects A and B to be (Conjugate)Transposed." );
    if( A.Width()  != C.Height() ||
        B.Height() != C.Width()  ||
        A.Height() != B.Width()    )
    {
        ostringstream msg;
        msg << "Nonconformal GemmTTB: " << endl
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

    DistMatrix<T,MC,MR> CT(g),  C0(g),
                        CB(g),  C1(g),
                                C2(g);

    // Temporary distributions
    DistMatrix<T,MR,  Star> A1_MR_Star(g);
    DistMatrix<T,Star,MC  > D1_Star_MC(g);
    DistMatrix<T,MR,  MC  > D1_MR_MC(g);
    DistMatrix<T,MC,  MR  > D1(g);

    // Start the algorithm 
    blas::Scal( beta, C );
    LockedPartitionRight( A, AL, AR, 0 );
    PartitionDown
    ( C, CT,
         CB, 0 );
    while( AR.Width() > 0 )
    {
        LockedRepartitionRight
        ( AL, /**/     AR,
          A0, /**/ A1, A2 );
 
        RepartitionDown
        ( CT,  C0,
         /**/ /**/
               C1,
          CB,  C2 );

        A1_MR_Star.AlignWith( B );
        D1_Star_MC.AlignWith( B );
        D1_Star_MC.ResizeTo( C1.Height(), C1.Width() );
        D1.AlignWith( C1 );
        //--------------------------------------------------------------------//
        A1_MR_Star = A1; // A1[MR,*] <- A1[MC,MR]
 
        // D1[*,MC] := alpha (A1[MR,*])^T (B[MC,MR])^T
        //           = alpha (A1^T)[*,MR] (B^T)[MR,MC]
        blas::internal::LocalGemm
        ( orientationOfA, orientationOfB, 
          alpha, A1_MR_Star, B, (T)0, D1_Star_MC );

        // C1[MC,MR] += scattered & transposed D1[*,MC] summed over grid rows
        D1_MR_MC.SumScatterFrom( D1_Star_MC );
        D1 = D1_MR_MC; 
        blas::Axpy( (T)1, D1, C1 );
        //--------------------------------------------------------------------//
        A1_MR_Star.FreeAlignments();
        D1_Star_MC.FreeAlignments();
        D1.FreeAlignments();

        SlideLockedPartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );

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

// Transpose Transpose Gemm that avoids communicating the matrix C.
template<typename T>
void
elemental::blas::internal::GemmTTC
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::GemmTTC");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
    if( orientationOfA == Normal || orientationOfB == Normal )
        throw logic_error
        ( "GemmTTC expects A and B to be (Conjugate)Transposed." );
    if( A.Width()  != C.Height() ||
        B.Height() != C.Width()  ||
        A.Height() != B.Width()    )
    {
        ostringstream msg;
        msg << "Nonconformal GemmTTC: " << endl
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

    DistMatrix<T,MC,MR> BL(g), BR(g),
                        B0(g), B1(g), B2(g);

    // Temporary distributions
    DistMatrix<T,Star,MC> A1_Star_MC(g);
    DistMatrix<T,MR,Star> B1_MR_Star(g);
    
    // Start the algorithm    
    blas::Scal( beta, C );
    LockedPartitionDown
    ( A, AT,
         AB, 0 ); 
    LockedPartitionRight( B, BL, BR, 0 );
    while( AB.Height() > 0 )
    {
        LockedRepartitionDown
        ( AT,  A0,
         /**/ /**/
               A1,
          AB,  A2 );

        LockedRepartitionRight
        ( BL, /**/     BR,
          B0, /**/ B1, B2 );

        A1_Star_MC.AlignWith( C );
        B1_MR_Star.AlignWith( C );
        //--------------------------------------------------------------------//
        A1_Star_MC = A1; // A1[*,MC] <- A1[MC,MR]
        B1_MR_Star = B1; // B1[MR,*] <- B1[MC,MR]

        // C[MC,MR] += alpha (A1[*,MC])^T (B1[MR,*])^T
        //           = alpha (A1^T)[MC,*] (B1^T)[*,MR]
        blas::internal::LocalGemm
        ( orientationOfA, orientationOfB, 
          alpha, A1_Star_MC, B1_MR_Star, (T)1, C );
        //--------------------------------------------------------------------//
        A1_Star_MC.FreeAlignments();
        B1_MR_Star.FreeAlignments();

        SlideLockedPartitionDown
        ( AT,  A0,
               A1,
         /**/ /**/
          AB,  A2 );

        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::GemmTT
( Orientation orientationOfA, 
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,MC,MR>& A, 
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::GemmTTA
( Orientation orientationOfA, 
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,MC,MR>& A, 
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::GemmTTB
( Orientation orientationOfA, 
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,MC,MR>& A, 
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::GemmTTC
( Orientation orientationOfA, 
  Orientation orientationOfB,
  float alpha, const DistMatrix<float,MC,MR>& A, 
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::GemmTT
( Orientation orientationOfA, 
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,MC,MR>& A, 
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

template void elemental::blas::internal::GemmTTA
( Orientation orientationOfA, 
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,MC,MR>& A, 
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

template void elemental::blas::internal::GemmTTB
( Orientation orientationOfA, 
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,MC,MR>& A, 
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

template void elemental::blas::internal::GemmTTC
( Orientation orientationOfA, 
  Orientation orientationOfB,
  double alpha, const DistMatrix<double,MC,MR>& A, 
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::GemmTT
( Orientation orientationOfA, 
  Orientation orientationOfB,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A, 
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::GemmTTA
( Orientation orientationOfA, 
  Orientation orientationOfB,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A, 
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::GemmTTB
( Orientation orientationOfA, 
  Orientation orientationOfB,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A, 
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::GemmTTC
( Orientation orientationOfA, 
  Orientation orientationOfB,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A, 
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::GemmTT
( Orientation orientationOfA, 
  Orientation orientationOfB,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A, 
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void elemental::blas::internal::GemmTTA
( Orientation orientationOfA, 
  Orientation orientationOfB,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A, 
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void elemental::blas::internal::GemmTTB
( Orientation orientationOfA, 
  Orientation orientationOfB,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A, 
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void elemental::blas::internal::GemmTTC
( Orientation orientationOfA, 
  Orientation orientationOfB,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A, 
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

