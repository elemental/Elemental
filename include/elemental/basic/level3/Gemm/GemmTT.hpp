/*
   Copyright (c) 2009-2011, Jack Poulson
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

template<typename T>
inline void
elemental::basic::internal::GemmTT
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::GemmTT");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
    if( orientationOfA == NORMAL || orientationOfB == NORMAL )
        throw std::logic_error("GemmTT expects A and B to be transposed");
#endif
    const int m = C.Height();
    const int n = C.Width();
    const int k = A.Height();
    const float weightTowardsC = 2.0;

    if( m <= n && weightTowardsC*m <= k )
    {
        basic::internal::GemmTTB
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
    else if( n <= m && weightTowardsC*n <= k )
    {
        basic::internal::GemmTTA
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
    else
    {
        basic::internal::GemmTTC
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Transpose Transpose Gemm that avoids communicating the matrix A.
template<typename T>
inline void
elemental::basic::internal::GemmTTA
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::GemmTTA");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
    if( orientationOfA == NORMAL || orientationOfB == NORMAL )
        throw std::logic_error
        ("GemmTTA expects A and B to be (Conjugate)Transposed");
    if( A.Width()  != C.Height() ||
        B.Height() != C.Width()  ||
        A.Height() != B.Width()    )
    {
        std::ostringstream msg;
        msg << "Nonconformal GemmTTA: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  B ~ " << B.Height() << " x " << B.Width() << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T,MC,MR> BT(g),  B0(g),
                        BB(g),  B1(g),
                                B2(g);

    DistMatrix<T,MC,MR> CL(g), CR(g),
                        C0(g), C1(g), C2(g);

    // Temporary distributions
    DistMatrix<T,STAR,MC  > B1_STAR_MC(g);
    DistMatrix<T,MR,  STAR> D1_MR_STAR(g);
    DistMatrix<T,MR,  MC  > D1_MR_MC(g);
    DistMatrix<T,MC,  MR  > D1(g);

    // Start the algorithm
    basic::Scal( beta, C );
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

        B1_STAR_MC.AlignWith( A ); 
        D1_MR_STAR.AlignWith( A );  
        D1_MR_STAR.ResizeTo( C1.Height(), C1.Width() );
        D1.AlignWith( C1 );  
        //--------------------------------------------------------------------//
        B1_STAR_MC = B1; // B1[*,MC] <- B1[MC,MR]

        // D1[MR,*] := alpha (A[MC,MR])^T (B1[*,MC])^T
        //           = alpha (A^T)[MR,MC] (B1^T)[MC,*]
        basic::internal::LocalGemm
        ( orientationOfA, orientationOfB, 
          alpha, A, B1_STAR_MC, (T)0, D1_MR_STAR );

        // C1[MC,MR] += scattered & transposed D1[MR,*] summed over grid cols
        D1_MR_MC.SumScatterFrom( D1_MR_STAR );
        D1 = D1_MR_MC; 
        basic::Axpy( (T)1, D1, C1 );
        //--------------------------------------------------------------------//
        B1_STAR_MC.FreeAlignments();
        D1_MR_STAR.FreeAlignments();
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
inline void
elemental::basic::internal::GemmTTB
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::GemmTTB");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
    if( orientationOfA == NORMAL || orientationOfB == NORMAL )
        throw std::logic_error
        ("GemmTTB expects A and B to be (Conjugate)Transposed");
    if( A.Width()  != C.Height() ||
        B.Height() != C.Width()  ||
        A.Height() != B.Width()    )
    {
        std::ostringstream msg;
        msg << "Nonconformal GemmTTB: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  B ~ " << B.Height() << " x " << B.Width() << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T,MC,MR> AL(g), AR(g),
                        A0(g), A1(g), A2(g);

    DistMatrix<T,MC,MR> CT(g),  C0(g),
                        CB(g),  C1(g),
                                C2(g);

    // Temporary distributions
    DistMatrix<T,MR,  STAR> A1_MR_STAR(g);
    DistMatrix<T,STAR,MC  > D1_STAR_MC(g);
    DistMatrix<T,MR,  MC  > D1_MR_MC(g);
    DistMatrix<T,MC,  MR  > D1(g);

    // Start the algorithm 
    basic::Scal( beta, C );
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

        A1_MR_STAR.AlignWith( B );
        D1_STAR_MC.AlignWith( B );
        D1_STAR_MC.ResizeTo( C1.Height(), C1.Width() );
        D1.AlignWith( C1 );
        //--------------------------------------------------------------------//
        A1_MR_STAR = A1; // A1[MR,*] <- A1[MC,MR]
 
        // D1[*,MC] := alpha (A1[MR,*])^T (B[MC,MR])^T
        //           = alpha (A1^T)[*,MR] (B^T)[MR,MC]
        basic::internal::LocalGemm
        ( orientationOfA, orientationOfB, 
          alpha, A1_MR_STAR, B, (T)0, D1_STAR_MC );

        // C1[MC,MR] += scattered & transposed D1[*,MC] summed over grid rows
        D1_MR_MC.SumScatterFrom( D1_STAR_MC );
        D1 = D1_MR_MC; 
        basic::Axpy( (T)1, D1, C1 );
        //--------------------------------------------------------------------//
        A1_MR_STAR.FreeAlignments();
        D1_STAR_MC.FreeAlignments();
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
inline void
elemental::basic::internal::GemmTTC
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("basic::internal::GemmTTC");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
    if( orientationOfA == NORMAL || orientationOfB == NORMAL )
        throw std::logic_error
        ("GemmTTC expects A and B to be (Conjugate)Transposed");
    if( A.Width()  != C.Height() ||
        B.Height() != C.Width()  ||
        A.Height() != B.Width()    )
    {
        std::ostringstream msg;
        msg << "Nonconformal GemmTTC: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  B ~ " << B.Height() << " x " << B.Width() << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str() );
    }
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T,MC,MR> AT(g),  A0(g),
                        AB(g),  A1(g),
                                A2(g);

    DistMatrix<T,MC,MR> BL(g), BR(g),
                        B0(g), B1(g), B2(g);

    // Temporary distributions
    DistMatrix<T,STAR,MC> A1_STAR_MC(g);
    DistMatrix<T,MR,STAR> B1_MR_STAR(g);
    
    // Start the algorithm    
    basic::Scal( beta, C );
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

        A1_STAR_MC.AlignWith( C );
        B1_MR_STAR.AlignWith( C );
        //--------------------------------------------------------------------//
        A1_STAR_MC = A1; // A1[*,MC] <- A1[MC,MR]
        B1_MR_STAR = B1; // B1[MR,*] <- B1[MC,MR]

        // C[MC,MR] += alpha (A1[*,MC])^T (B1[MR,*])^T
        //           = alpha (A1^T)[MC,*] (B1^T)[*,MR]
        basic::internal::LocalGemm
        ( orientationOfA, orientationOfB, 
          alpha, A1_STAR_MC, B1_MR_STAR, (T)1, C );
        //--------------------------------------------------------------------//
        A1_STAR_MC.FreeAlignments();
        B1_MR_STAR.FreeAlignments();

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
