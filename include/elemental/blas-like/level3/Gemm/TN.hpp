/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_GEMM_TN_HPP
#define BLAS_GEMM_TN_HPP

namespace elem {
namespace internal {

// Transpose Normal Gemm that avoids communicating the matrix A.
template<typename T> 
inline void
GemmTNA
( Orientation orientationOfA,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmTNA");    
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
    if( orientationOfA == NORMAL )
        throw std::logic_error("GemmTNA assumes A is (Conjugate)Transposed");
    if( A.Width()  != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Height() != B.Height()   )
    {
        std::ostringstream msg;
        msg << "Nonconformal GemmTNA: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  B ~ " << B.Height() << " x " << B.Width() << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T> BL(g), BR(g),
                  B0(g), B1(g), B2(g);
    DistMatrix<T> CL(g), CR(g),
                  C0(g), C1(g), C2(g);

    // Temporary distributions
    DistMatrix<T,MC,STAR> B1_MC_STAR(g);
    DistMatrix<T,MR,STAR> D1_MR_STAR(g);
    DistMatrix<T,MR,MC  > D1_MR_MC(g);
    DistMatrix<T> D1(g);

    B1_MC_STAR.AlignWith( A );
    D1_MR_STAR.AlignWith( A );

    // Start the algorithm
    Scale( beta, C );
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

        D1.AlignWith( C1 );
        Zeros( C1.Height(), C1.Width(), D1_MR_STAR );
        //--------------------------------------------------------------------//
        B1_MC_STAR = B1; // B1[MC,*] <- B1[MC,MR]

        // D1[MR,*] := alpha (A1[MC,MR])^T B1[MC,*]
        //           = alpha (A1^T)[MR,MC] B1[MC,*]
        LocalGemm
        ( orientationOfA, NORMAL, alpha, A, B1_MC_STAR, T(0), D1_MR_STAR );

        // C1[MC,MR] += scattered & transposed D1[MR,*] summed over grid cols
        D1_MR_MC.SumScatterFrom( D1_MR_STAR );
        D1 = D1_MR_MC; 
        Axpy( T(1), D1, C1 );
        //--------------------------------------------------------------------//
        D1.FreeAlignments();

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

// Transpose Normal Gemm that avoids communicating the matrix B.
template<typename T> 
inline void
GemmTNB
( Orientation orientationOfA,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmTNB");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
    if( orientationOfA == NORMAL )
        throw std::logic_error("GemmTNB assumes A is (Conjugate)Transposed");
    if( A.Width()  != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Height() != B.Height()   )
    {
        std::ostringstream msg;
        msg << "Nonconformal GemmTNB: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  B ~ " << B.Height() << " x " << B.Width() << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str() );
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
    DistMatrix<T,MC,STAR> A1_MC_STAR(g);
    DistMatrix<T,MR,STAR> D1AdjOrTrans_MR_STAR(g);

    A1_MC_STAR.AlignWith( B );
    D1AdjOrTrans_MR_STAR.AlignWith( B );

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

        Zeros( C1.Width(), C1.Height(), D1AdjOrTrans_MR_STAR );
        //--------------------------------------------------------------------//
        A1_MC_STAR = A1; // A1[MC,*] <- A1[MC,MR]

        // D1[*,MR] := alpha (A1[MC,*])^[T/H] B[MC,MR]
        //           = alpha (A1^[T/H])[*,MC] B[MC,MR]
        if( orientationOfA == ADJOINT )
        {
            LocalGemm
            ( ADJOINT, NORMAL, 
              Conj(alpha), B, A1_MC_STAR, T(0), D1AdjOrTrans_MR_STAR );
            C1.AdjointSumScatterUpdate( T(1), D1AdjOrTrans_MR_STAR );
        }
        else
        {
            LocalGemm
            ( TRANSPOSE, NORMAL,
              alpha, B, A1_MC_STAR, T(0), D1AdjOrTrans_MR_STAR );
            C1.TransposeSumScatterUpdate( T(1), D1AdjOrTrans_MR_STAR );
        }
        //--------------------------------------------------------------------//

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

// Transpose Normal Gemm that avoids communicating the matrix C.
template<typename T> 
inline void
GemmTNC
( Orientation orientationOfA,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmTNC");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
    if( orientationOfA == NORMAL )
        throw std::logic_error("GemmTNC assumes A is (Conjugate)Transposed");
    if( A.Width()  != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Height() != B.Height()   )
    {
        std::ostringstream msg;
        msg << "Nonconformal GemmTNC: \n"
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
    DistMatrix<T> BT(g),  B0(g),
                  BB(g),  B1(g),
                          B2(g);

    // Temporary distributions
    DistMatrix<T,STAR,MC> A1_STAR_MC(g);
    DistMatrix<T,MR,STAR> B1Trans_MR_STAR(g);

    A1_STAR_MC.AlignWith( C );
    B1Trans_MR_STAR.AlignWith( C );

    // Start the algorithm
    Scale( beta, C );
    LockedPartitionDown
    ( A, AT,
         AB, 0 );
    LockedPartitionDown
    ( B, BT,
         BB, 0 );
    while( AB.Height() > 0 )
    {
        LockedRepartitionDown
        ( AT,  A0,
         /**/ /**/
               A1,
          AB,  A2 );

        LockedRepartitionDown
        ( BT,  B0,
         /**/ /**/
               B1,
          BB,  B2 );

        //--------------------------------------------------------------------//
        A1_STAR_MC = A1; 
        B1Trans_MR_STAR.TransposeFrom( B1 );

        // C[MC,MR] += alpha (A1[*,MC])^T B1[*,MR]
        //           = alpha (A1^T)[MC,*] B1[*,MR]
        LocalGemm
        ( orientationOfA, TRANSPOSE, 
          alpha, A1_STAR_MC, B1Trans_MR_STAR, T(1), C );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDown
        ( AT,  A0,
               A1,
         /**/ /**/
          AB,  A2 );

        SlideLockedPartitionDown
        ( BT,  B0,
               B1,
         /**/ /**/
          BB,  B2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
GemmTN
( Orientation orientationOfA,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmTN");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
    if( orientationOfA == NORMAL )
        throw std::logic_error("GemmTN assumes A is (Conjugate)Transposed");
#endif
    const int m = C.Height();
    const int n = C.Width();
    const int k = A.Height();
    const float weightTowardsC = 2.0;

    if( m <= n && weightTowardsC*m <= k )
    {
        GemmTNB( orientationOfA, alpha, A, B, beta, C );
    }
    else if( n <= m && weightTowardsC*n <= k )
    {
        GemmTNA( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        GemmTNC( orientationOfA, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem

#endif // ifndef BLAS_GEMM_TN_HPP
