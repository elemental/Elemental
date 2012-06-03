/*
   Copyright (c) 2009-2012, Jack Poulson
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

namespace elem {
namespace internal {

template<typename T>
inline void
GemmNT
( Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmNT");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
    if( orientationOfB == NORMAL )
        throw std::logic_error
        ("GemmNT requires that B be (Conjugate)Transposed");
#endif
    const int m = C.Height();
    const int n = C.Width();
    const int k = A.Width();
    const float weightTowardsC = 2.0;

    if( m <= n && weightTowardsC*m <= k )
    {
        GemmNTB( orientationOfB, alpha, A, B, beta, C );
    }
    else if( n <= m && weightTowardsC*n <= k )
    {
        GemmNTA( orientationOfB, alpha, A, B, beta, C );
    }
    else
    {
        GemmNTC( orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Normal Transpose Gemm that avoids communicating the matrix A.
template<typename T>
inline void
GemmNTA
( Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmNTA");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
    if( orientationOfB == NORMAL )
        throw std::logic_error
        ("GemmTNA requires that B be (Conjugate)Transposed");
    if( A.Height() != C.Height() ||
        B.Height() != C.Width()  ||
        A.Width()  != B.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal GemmNTA: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  B ~ " << B.Height() << " x " << B.Width() << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T> BT(g),  B0(g),
                  BB(g),  B1(g),
                          B2(g);

    DistMatrix<T> CL(g), CR(g),
                  C0(g), C1(g), C2(g);

    // Temporary distributions
    DistMatrix<T,STAR,MR> B1_STAR_MR(g);
    DistMatrix<T,MC,STAR> D1_MC_STAR(g);

    // Start the algorithm
    Scale( beta, C );
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

        B1_STAR_MR.AlignWith( A );
        D1_MC_STAR.AlignWith( A );
        D1_MC_STAR.ResizeTo( C1.Height(), C1.Width() );
        //--------------------------------------------------------------------//
        B1_STAR_MR = B1; // B1[*,MR] <- B1[MC,MR]

        // C1[MC,*] := alpha A[MC,MR] (B1[*,MR])^T
        //           = alpha A[MC,MR] (B1^T)[MR,*]
        LocalGemm
        ( NORMAL, orientationOfB, alpha, A, B1_STAR_MR, (T)0, D1_MC_STAR );

        // C1[MC,MR] += scattered result of D1[MC,*] summed over grid rows
        C1.SumScatterUpdate( (T)1, D1_MC_STAR );
        //--------------------------------------------------------------------//
        B1_STAR_MR.FreeAlignments();
        D1_MC_STAR.FreeAlignments();

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
inline void
GemmNTB
( Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmNTB");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
    if( orientationOfB == NORMAL )
        throw std::logic_error
        ("GemmNTB requires that B be (Conjugate)Transposed");
    if( A.Height() != C.Height() ||
        B.Height() != C.Width()  ||
        A.Width()  != B.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal GemmNTB: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  B ~ " << B.Height() << " x " << B.Width() << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T> AT(g),  A0(g),
                  AB(g),  A1(g),
                          A2(g);

    DistMatrix<T> CT(g),  C0(g),
                  CB(g),  C1(g),
                          C2(g);

    // Temporary distributions
    DistMatrix<T,STAR,MR> A1_STAR_MR(g);
    DistMatrix<T,STAR,MC> D1_STAR_MC(g);
    DistMatrix<T,MR,  MC> D1_MR_MC(g);
    DistMatrix<T> D1(g);

    // Start the algorithm
    Scale( beta, C );
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

        A1_STAR_MR.AlignWith( B );
        D1_STAR_MC.AlignWith( B );
        D1_STAR_MC.ResizeTo( C1.Height(), C1.Width() );
        D1.AlignWith( C1 );
        //--------------------------------------------------------------------//
        A1_STAR_MR = A1; // A1[*,MR] <- A1[MC,MR]

        // D1[*,MC] := alpha A1[*,MR] (B[MC,MR])^T
        //           = alpha A1[*,MR] (B^T)[MR,MC]
        LocalGemm
        ( NORMAL, orientationOfB, alpha, A1_STAR_MR, B, (T)0, D1_STAR_MC );

        // C1[MC,MR] += scattered & transposed D1[*,MC] summed over grid rows
        D1_MR_MC.SumScatterFrom( D1_STAR_MC );
        D1 = D1_MR_MC; 
        Axpy( (T)1, D1, C1 );
        //--------------------------------------------------------------------//
        A1_STAR_MR.FreeAlignments();
        D1_STAR_MC.FreeAlignments();
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
inline void
GemmNTC
( Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmNTC");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
    if( orientationOfB == NORMAL )
        throw std::logic_error
        ("GemmNTC requires that B be (Conjugate)Transposed");
    if( A.Height() != C.Height() ||
        B.Height() != C.Width()  ||
        A.Width()  != B.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal GemmNTC: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  B ~ " << B.Height() << " x " << B.Width() << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T> AL(g), AR(g),
                  A0(g), A1(g), A2(g);

    DistMatrix<T> BL(g), BR(g),
                  B0(g), B1(g), B2(g);

    // Temporary distributions
    DistMatrix<T,MC,STAR> A1_MC_STAR(g);
    DistMatrix<T,MR,STAR> B1_MR_STAR(g);

    // Start the algorithm
    Scale( beta, C );
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

        A1_MC_STAR.AlignWith( C );
        B1_MR_STAR.AlignWith( C );
        //--------------------------------------------------------------------//
        A1_MC_STAR = A1; // A1[MC,*] <- A1[MC,MR]
        B1_MR_STAR = B1; // B1[MR,*] <- B1[MC,MR]

        // C[MC,MR] += alpha A1[MC,*] (B1[MR,*])^T
        //           = alpha A1[MC,*] (B1^T)[*,MR]
        LocalGemm
        ( NORMAL, orientationOfB, alpha, A1_MC_STAR, B1_MR_STAR, (T)1, C );
        //--------------------------------------------------------------------//
        A1_MC_STAR.FreeAlignments();
        B1_MR_STAR.FreeAlignments();
 
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

} // namespace internal
} // namespace elem
