/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

// Transpose Transpose Gemm that avoids communicating the matrix A.
template<typename T>
inline void
GemmTTA
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmTTA");
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
    DistMatrix<T> BT(g),  B0(g),
                  BB(g),  B1(g),
                          B2(g);
    DistMatrix<T> CL(g), CR(g),
                  C0(g), C1(g), C2(g);

    // Temporary distributions
    DistMatrix<T,STAR,MC  > B1_STAR_MC(g);
    DistMatrix<T,MR,  STAR> D1_MR_STAR(g);
    DistMatrix<T,MR,  MC  > D1_MR_MC(g);
    DistMatrix<T> D1(g);

    B1_STAR_MC.AlignWith( A ); 
    D1_MR_STAR.AlignWith( A );  

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

        D1.AlignWith( C1 );  
        Zeros( C1.Height(), C1.Width(), D1_MR_STAR );
        //--------------------------------------------------------------------//
        B1_STAR_MC = B1; // B1[*,MC] <- B1[MC,MR]

        // D1[MR,*] := alpha (A[MC,MR])^T (B1[*,MC])^T
        //           = alpha (A^T)[MR,MC] (B1^T)[MC,*]
        LocalGemm
        ( orientationOfA, orientationOfB, 
          alpha, A, B1_STAR_MC, T(0), D1_MR_STAR );

        // C1[MC,MR] += scattered & transposed D1[MR,*] summed over grid cols
        D1_MR_MC.SumScatterFrom( D1_MR_STAR );
        D1 = D1_MR_MC; 
        Axpy( T(1), D1, C1 );
        //--------------------------------------------------------------------//
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
GemmTTB
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmTTB");
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
    DistMatrix<T> AL(g), AR(g),
                  A0(g), A1(g), A2(g);
    DistMatrix<T> CT(g),  C0(g),
                  CB(g),  C1(g),
                          C2(g);

    // Temporary distributions
    DistMatrix<T,VR,  STAR> A1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > A1AdjOrTrans_STAR_MR(g);
    DistMatrix<T,STAR,MC  > D1_STAR_MC(g);
    DistMatrix<T,MR,  MC  > D1_MR_MC(g);
    DistMatrix<T> D1(g);

    A1_VR_STAR.AlignWith( B );
    A1AdjOrTrans_STAR_MR.AlignWith( B );
    D1_STAR_MC.AlignWith( B );

    // Start the algorithm 
    Scale( beta, C );
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

        D1.AlignWith( C1 );
        Zeros( C1.Height(), C1.Width(), D1_STAR_MC );
        //--------------------------------------------------------------------//
        A1_VR_STAR = A1;
        if( orientationOfA == ADJOINT )
            A1AdjOrTrans_STAR_MR.AdjointFrom( A1_VR_STAR );
        else
            A1AdjOrTrans_STAR_MR.TransposeFrom( A1_VR_STAR );
 
        // D1[*,MC] := alpha (A1[MR,*])^[T/H] (B[MC,MR])^[T/H]
        //           = alpha (A1^[T/H])[*,MR] (B^[T/H])[MR,MC]
        LocalGemm
        ( NORMAL, orientationOfB, 
          alpha, A1AdjOrTrans_STAR_MR, B, T(0), D1_STAR_MC );

        // C1[MC,MR] += scattered & transposed D1[*,MC] summed over grid rows
        D1_MR_MC.SumScatterFrom( D1_STAR_MC );
        D1 = D1_MR_MC; 
        Axpy( T(1), D1, C1 );
        //--------------------------------------------------------------------//
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
GemmTTC
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmTTC");
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
    DistMatrix<T> AT(g),  A0(g),
                  AB(g),  A1(g),
                          A2(g);
    DistMatrix<T> BL(g), BR(g),
                  B0(g), B1(g), B2(g);

    // Temporary distributions
    DistMatrix<T,STAR,MC  > A1_STAR_MC(g);
    DistMatrix<T,VR,  STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > B1AdjOrTrans_STAR_MR(g);

    A1_STAR_MC.AlignWith( C );
    B1_VR_STAR.AlignWith( C );
    B1AdjOrTrans_STAR_MR.AlignWith( C );
    
    // Start the algorithm    
    Scale( beta, C );
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

        //--------------------------------------------------------------------//
        A1_STAR_MC = A1; 
        B1_VR_STAR = B1;
        if( orientationOfB == ADJOINT )
            B1AdjOrTrans_STAR_MR.AdjointFrom( B1_VR_STAR );
        else
            B1AdjOrTrans_STAR_MR.TransposeFrom( B1_VR_STAR );

        // C[MC,MR] += alpha (A1[*,MC])^[T/H] (B1[MR,*])^[T/H]
        //           = alpha (A1^[T/H])[MC,*] (B1^[T/H])[*,MR]
        LocalGemm
        ( orientationOfA, NORMAL, 
          alpha, A1_STAR_MC, B1AdjOrTrans_STAR_MR, T(1), C );
        //--------------------------------------------------------------------//

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

template<typename T>
inline void
GemmTT
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmTT");
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
        GemmTTB( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
    else if( n <= m && weightTowardsC*n <= k )
    {
        GemmTTA( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
    else
    {
        GemmTTC( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
