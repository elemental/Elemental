/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace gemm {

// Normal Transpose Gemm that avoids communicating the matrix A
template<typename T>
inline void
SUMMA_NTA
( Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("gemm::SUMMA_NTA");
        if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
            LogicError("{A,B,C} must have the same grid");
        if( orientationOfB == NORMAL )
            LogicError("B must be (Conjugate)Transposed");
        if( A.Height() != C.Height() || B.Height() != C.Width() ||
            A.Width() != B.Width() )
            LogicError
            ("Nonconformal matrices:\n",
             DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
    )
    const Grid& g = A.Grid();
    const bool conjugate = ( orientationOfB == ADJOINT );

    // Matrix views
    DistMatrix<T> BT(g),  B0(g),
                  BB(g),  B1(g),
                          B2(g);
    DistMatrix<T> CL(g), CR(g),
                  C0(g), C1(g), C2(g);

    // Temporary distributions
    DistMatrix<T,MR,STAR> B1Trans_MR_STAR(g);
    DistMatrix<T,MC,STAR> D1_MC_STAR(g);

    B1Trans_MR_STAR.AlignWith( A );
    D1_MC_STAR.AlignWith( A );

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

        //--------------------------------------------------------------------//
        B1.TransposeColAllGather( B1Trans_MR_STAR, conjugate );

        // C1[MC,*] := alpha A[MC,MR] (B1^[T/H])[MR,*]
        LocalGemm( NORMAL, NORMAL, alpha, A, B1Trans_MR_STAR, D1_MC_STAR );

        // C1[MC,MR] += scattered result of D1[MC,*] summed over grid rows
        C1.RowSumScatterUpdate( T(1), D1_MC_STAR );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDown
        ( BT,  B0,
               B1,
         /**/ /**/
          BB,  B2 );

        SlidePartitionRight
        ( CL,     /**/ CR,
          C0, C1, /**/ C2 );
    }
}

// Normal Transpose Gemm that avoids communicating the matrix B
template<typename T>
inline void
SUMMA_NTB
( Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("gemm::SUMMA_NTB");
        if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
            LogicError("{A,B,C} must have the same grid");
        if( orientationOfB == NORMAL )
            LogicError("B must be (Conjugate)Transposed");
        if( A.Height() != C.Height() ||
            B.Height() != C.Width() ||
            A.Width() != B.Width() )
            LogicError
            ("Nonconformal matrices:\n",
             DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
    )
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T> AT(g),  A0(g),
                  AB(g),  A1(g),
                          A2(g);
    DistMatrix<T> CT(g),  C0(g),
                  CB(g),  C1(g),
                          C2(g);

    // Temporary distributions
    DistMatrix<T,MR,STAR> A1Trans_MR_STAR(g);
    DistMatrix<T,STAR,MC> D1_STAR_MC(g);
    DistMatrix<T,MR,MC> D1_MR_MC(g);
    DistMatrix<T> D1(g);

    A1Trans_MR_STAR.AlignWith( B );
    D1_STAR_MC.AlignWith( B );

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

        D1.AlignWith( C1 );
        //--------------------------------------------------------------------//
        A1.TransposeColAllGather( A1Trans_MR_STAR );

        // D1[*,MC] := alpha A1[*,MR] (B[MC,MR])^T
        //           = alpha (A1^T)[MR,*] (B^T)[MR,MC]
        LocalGemm
        ( TRANSPOSE, orientationOfB, alpha, A1Trans_MR_STAR, B, D1_STAR_MC );

        // C1[MC,MR] += scattered & transposed D1[*,MC] summed over grid rows
        D1_MR_MC.ColSumScatterFrom( D1_STAR_MC );
        D1 = D1_MR_MC; 
        Axpy( T(1), D1, C1 );
        //--------------------------------------------------------------------//

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

// Normal Transpose Gemm that avoids communicating the matrix C
template<typename T>
inline void
SUMMA_NTC
( Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("gemm::SUMMA_NTC");
        if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
            LogicError("{A,B,C} must have the same grid");
        if( orientationOfB == NORMAL )
            LogicError("B must be (Conjugate)Transposed");
        if( A.Height() != C.Height() ||
            B.Height() != C.Width() ||
            A.Width() != B.Width() )
            LogicError
            ("Nonconformal matrices:\n",
             DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
    )
    const Grid& g = A.Grid();
    const bool conjugate = ( orientationOfB == ADJOINT );

    // Matrix views
    DistMatrix<T> AL(g), AR(g),
                  A0(g), A1(g), A2(g);
    DistMatrix<T> BL(g), BR(g),
                  B0(g), B1(g), B2(g);

    // Temporary distributions
    DistMatrix<T,MC,STAR> A1_MC_STAR(g);
    DistMatrix<T,VR,STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR> B1Trans_STAR_MR(g);

    A1_MC_STAR.AlignWith( C );
    B1_VR_STAR.AlignWith( C );
    B1Trans_STAR_MR.AlignWith( C );

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

        //--------------------------------------------------------------------//
        A1_MC_STAR = A1; // A1[MC,*] <- A1[MC,MR]
        B1_VR_STAR = B1;
        B1_VR_STAR.TransposePartialColAllGather( B1Trans_STAR_MR, conjugate );

        // C[MC,MR] += alpha A1[MC,*] (B1[MR,*])^T
        LocalGemm
        ( NORMAL, NORMAL, alpha, A1_MC_STAR, B1Trans_STAR_MR, T(1), C );
        //--------------------------------------------------------------------//
 
        SlideLockedPartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );

        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );
    }
}

template<typename T>
inline void
SUMMA_NT
( Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta, DistMatrix<T>& C, GemmAlgorithm alg=GEMM_DEFAULT )
{
    DEBUG_ONLY(CallStackEntry cse("gemm::SUMMA_NT"))
    const Int m = C.Height();
    const Int n = C.Width();
    const Int k = A.Width();
    const double weightTowardsC = 2.;

    switch( alg )
    {
    case GEMM_DEFAULT:
        if( m <= n && weightTowardsC*m <= k )
            SUMMA_NTB( orientationOfB, alpha, A, B, beta, C );
        else if( n <= m && weightTowardsC*n <= k )
            SUMMA_NTA( orientationOfB, alpha, A, B, beta, C );
        else
            SUMMA_NTC( orientationOfB, alpha, A, B, beta, C );
        break;
    case GEMM_SUMMA_A: SUMMA_NTA( orientationOfB, alpha, A, B, beta, C ); break;
    case GEMM_SUMMA_B: SUMMA_NTB( orientationOfB, alpha, A, B, beta, C ); break;
    case GEMM_SUMMA_C: SUMMA_NTC( orientationOfB, alpha, A, B, beta, C ); break;
    default: LogicError("Unsupported Gemm option");
    }
}

} // namespace gemm
} // namespace El
