/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace gemm {

// Transpose Normal Gemm that avoids communicating the matrix A
template<typename T> 
inline void
SUMMA_TNA
( Orientation orientationOfA,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("gemm::SUMMA_TNA");    
        if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
            LogicError("{A,B,C} must have the same grid");
        if( orientationOfA == NORMAL )
            LogicError("A must be (Conjugate)Transposed");
        if( A.Width() != C.Height() || B.Width() != C.Width() ||
            A.Height() != B.Height() )
            LogicError
            ("Nonconformal matrices:\n",
             DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
    )
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
        //--------------------------------------------------------------------//
        B1_MC_STAR = B1; // B1[MC,*] <- B1[MC,MR]

        // D1[MR,*] := alpha (A1[MC,MR])^T B1[MC,*]
        //           = alpha (A1^T)[MR,MC] B1[MC,*]
        LocalGemm( orientationOfA, NORMAL, alpha, A, B1_MC_STAR, D1_MR_STAR );

        // C1[MC,MR] += scattered & transposed D1[MR,*] summed over grid cols
        D1_MR_MC.RowSumScatterFrom( D1_MR_STAR );
        D1 = D1_MR_MC; 
        Axpy( T(1), D1, C1 );
        //--------------------------------------------------------------------//

        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );

        SlidePartitionRight
        ( CL,     /**/ CR,
          C0, C1, /**/ C2 );
    }
}

// Transpose Normal Gemm that avoids communicating the matrix B
template<typename T> 
inline void
SUMMA_TNB
( Orientation orientationOfA,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("gemm::SUMMA_TNB");
        if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
            LogicError("{A,B,C} must have the same grid");
        if( orientationOfA == NORMAL )
            LogicError("A must be (Conjugate)Transposed");
        if( A.Width() != C.Height() || B.Width() != C.Width() ||
            A.Height() != B.Height() )
            LogicError
            ("Nonconformal matrices:\n",
             DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
    )
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T> AL(g), AR(g),
                  A0(g), A1(g), A2(g);
    DistMatrix<T> CT(g),  C0(g),
                  CB(g),  C1(g),
                          C2(g);

    // Temporary distributions
    DistMatrix<T,MC,STAR> A1_MC_STAR(g);
    DistMatrix<T,MR,STAR> D1Trans_MR_STAR(g);

    A1_MC_STAR.AlignWith( B );
    D1Trans_MR_STAR.AlignWith( B );

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

        //--------------------------------------------------------------------//
        A1_MC_STAR = A1; // A1[MC,*] <- A1[MC,MR]

        // D1[*,MR] := alpha (A1[MC,*])^[T/H] B[MC,MR]
        //           = alpha (A1^[T/H])[*,MC] B[MC,MR]
        if( orientationOfA == ADJOINT )
        {
            LocalGemm
            ( ADJOINT, NORMAL, 
              Conj(alpha), B, A1_MC_STAR, D1Trans_MR_STAR );
            C1.AdjointColSumScatterUpdate( T(1), D1Trans_MR_STAR );
        }
        else
        {
            LocalGemm
            ( TRANSPOSE, NORMAL, alpha, B, A1_MC_STAR, D1Trans_MR_STAR );
            C1.TransposeColSumScatterUpdate( T(1), D1Trans_MR_STAR );
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
}

// Transpose Normal Gemm that avoids communicating the matrix C
template<typename T> 
inline void
SUMMA_TNC
( Orientation orientationOfA,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("gemm::SUMMA_TNC");
        if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
            LogicError("{A,B,C} must have the same grid");
        if( orientationOfA == NORMAL )
            LogicError("A must be (Conjugate)Transposed");
        if( A.Width() != C.Height() || B.Width() != C.Width() ||
            A.Height() != B.Height() )
            LogicError
            ("Nonconformal matrices:\n",
             DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
    )
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
        B1.TransposeColAllGather( B1Trans_MR_STAR );

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
}

template<typename T>
inline void
SUMMA_TN
( Orientation orientationOfA,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C, GemmAlgorithm alg=GEMM_DEFAULT )
{
    DEBUG_ONLY(CallStackEntry cse("gemm::SUMMA_TN"))
    const Int m = C.Height();
    const Int n = C.Width();
    const Int k = A.Height();
    const double weightTowardsC = 2.;

    switch( alg )
    {
    case GEMM_DEFAULT:
        if( m <= n && weightTowardsC*m <= k )
            SUMMA_TNB( orientationOfA, alpha, A, B, beta, C );
        else if( n <= m && weightTowardsC*n <= k )
            SUMMA_TNA( orientationOfA, alpha, A, B, beta, C );
        else
            SUMMA_TNC( orientationOfA, alpha, A, B, beta, C );
        break;
    case GEMM_SUMMA_A: SUMMA_TNA( orientationOfA, alpha, A, B, beta, C ); break;
    case GEMM_SUMMA_B: SUMMA_TNB( orientationOfA, alpha, A, B, beta, C ); break;
    case GEMM_SUMMA_C: SUMMA_TNC( orientationOfA, alpha, A, B, beta, C ); break;
    default: LogicError("Unsupported Gemm option");
    }
}

} // namespace gemm
} // namespace El
