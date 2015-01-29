/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace gemm {

// Transpose Transpose Gemm that avoids communicating the matrix A
template<typename T>
inline void
SUMMA_TTA
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const AbstractDistMatrix<T>& APre, const AbstractDistMatrix<T>& BPre,
  T beta,        AbstractDistMatrix<T>& CPre )
{
    DEBUG_ONLY(
        CallStackEntry cse("gemm::SUMMA_TTA");
        AssertSameGrids( APre, BPre, CPre );
        if( orientationOfA == NORMAL || orientationOfB == NORMAL )
            LogicError("A and B must be (Conjugate)Transposed");
        if( APre.Width() != CPre.Height() || BPre.Height() != CPre.Width() ||
            APre.Height() != BPre.Width() )
            LogicError
            ("Nonconformal matrices:\n",
             DimsString(APre,"A"),"\n",DimsString(BPre,"B"),"\n",
             DimsString(CPre,"C"));
    )
    const Int m = CPre.Height();
    const Int n = CPre.Width();
    const Int sumDim = APre.Height();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();

    auto APtr = ReadProxy<T,MC,MR>( &APre );      auto& A = *APtr;
    auto BPtr = ReadProxy<T,MC,MR>( &BPre );      auto& B = *BPtr;
    auto CPtr = ReadWriteProxy<T,MC,MR>( &CPre ); auto& C = *CPtr;

    // Temporary distributions
    DistMatrix<T,STAR,MC  > B1_STAR_MC(g);
    DistMatrix<T,MR,  STAR> D1_MR_STAR(g);
    DistMatrix<T,MR,  MC  > D1_MR_MC(g);

    B1_STAR_MC.AlignWith( A ); 
    D1_MR_STAR.AlignWith( A );  

    Scale( beta, C );
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        auto B1 = B( IR(k,k+nb), IR(0,sumDim) );
        auto C1 = C( IR(0,m),    IR(k,k+nb)   );

        B1_STAR_MC = B1; // B1[*,MC] <- B1[MC,MR]

        // D1[MR,*] := alpha (A[MC,MR])^T (B1[*,MC])^T
        //           = alpha (A^T)[MR,MC] (B1^T)[MC,*]
        LocalGemm
        ( orientationOfA, orientationOfB, alpha, A, B1_STAR_MC, D1_MR_STAR );

        // C1[MC,MR] += scattered & transposed D1[MR,*] summed over grid cols
        Conjugate( D1_MR_STAR, D1_MR_MC );
        Axpy( T(1), D1_MR_MC, C1 );
    }
}

// Transpose Transpose Gemm that avoids communicating the matrix B
template<typename T>
inline void
SUMMA_TTB
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const AbstractDistMatrix<T>& APre, const AbstractDistMatrix<T>& BPre,
  T beta,        AbstractDistMatrix<T>& CPre )
{
    DEBUG_ONLY(
        CallStackEntry cse("gemm::SUMMA_TTB");
        AssertSameGrids( APre, BPre, CPre );
        if( orientationOfA == NORMAL || orientationOfB == NORMAL )
            LogicError("A and B must be (Conjugate)Transposed");
        if( APre.Width() != CPre.Height() || BPre.Height() != CPre.Width() ||
            APre.Height() != BPre.Width() )
            LogicError
            ("Nonconformal matrices:\n",
             DimsString(APre,"A"),"\n",DimsString(BPre,"B"),"\n",
             DimsString(CPre,"C"));
    )
    const Int m = CPre.Height();
    const Int n = CPre.Width();
    const Int sumDim = APre.Height();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();
    const bool conjugateA = ( orientationOfA == ADJOINT ); 

    auto APtr = ReadProxy<T,MC,MR>( &APre );      auto& A = *APtr;
    auto BPtr = ReadProxy<T,MC,MR>( &BPre );      auto& B = *BPtr;
    auto CPtr = ReadWriteProxy<T,MC,MR>( &CPre ); auto& C = *CPtr;

    // Temporary distributions
    DistMatrix<T,VR,  STAR> A1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > A1Trans_STAR_MR(g);
    DistMatrix<T,STAR,MC  > D1_STAR_MC(g);
    DistMatrix<T,MR,  MC  > D1_MR_MC(g);

    A1_VR_STAR.AlignWith( B );
    A1Trans_STAR_MR.AlignWith( B );
    D1_STAR_MC.AlignWith( B );

    Scale( beta, C );
    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);
        auto A1 = A( IR(0,sumDim), IR(k,k+nb) );
        auto C1 = C( IR(k,k+nb),   IR(0,n)    );

        // D1[*,MC] := alpha (A1[MR,*])^[T/H] (B[MC,MR])^[T/H]
        //           = alpha (A1^[T/H])[*,MR] (B^[T/H])[MR,MC]
        A1_VR_STAR = A1;
        Transpose( A1_VR_STAR, A1Trans_STAR_MR, conjugateA );
        LocalGemm
        ( NORMAL, orientationOfB, alpha, A1Trans_STAR_MR, B, D1_STAR_MC );

        // C1[MC,MR] += scattered & transposed D1[*,MC] summed over grid rows
        Contract( D1_STAR_MC, D1_MR_MC );
        Axpy( T(1), D1_MR_MC, C1 );
    }
}

// Transpose Transpose Gemm that avoids communicating the matrix C
template<typename T>
inline void
SUMMA_TTC
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const AbstractDistMatrix<T>& APre, const AbstractDistMatrix<T>& BPre,
  T beta,        AbstractDistMatrix<T>& CPre )
{
    DEBUG_ONLY(
        CallStackEntry cse("gemm::SUMMA_TTC");
        AssertSameGrids( APre, BPre, CPre );
        if( orientationOfA == NORMAL || orientationOfB == NORMAL )
            LogicError("A and B must be (Conjugate)Transposed");
        if( APre.Width() != CPre.Height() || BPre.Height() != CPre.Width() ||
            APre.Height() != BPre.Width() )
            LogicError
            ("Nonconformal matrices:\n",
             DimsString(APre,"A"),"\n",DimsString(BPre,"B"),"\n",
             DimsString(CPre,"C"));
    )
    const Int m = CPre.Height();
    const Int n = CPre.Width();
    const Int sumDim = APre.Height();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();
    const bool conjugateB = ( orientationOfB == ADJOINT );

    auto APtr = ReadProxy<T,MC,MR>( &APre );      auto& A = *APtr;
    auto BPtr = ReadProxy<T,MC,MR>( &BPre );      auto& B = *BPtr;
    auto CPtr = ReadWriteProxy<T,MC,MR>( &CPre ); auto& C = *CPtr;

    // Temporary distributions
    DistMatrix<T,STAR,MC  > A1_STAR_MC(g);
    DistMatrix<T,VR,  STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > B1Trans_STAR_MR(g);

    A1_STAR_MC.AlignWith( C );
    B1_VR_STAR.AlignWith( C );
    B1Trans_STAR_MR.AlignWith( C );
    
    Scale( beta, C );
    for( Int k=0; k<sumDim; k+=bsize )
    {
        const Int nb = Min(bsize,sumDim-k);
        auto A1 = A( IR(k,k+nb), IR(0,m)    );
        auto B1 = B( IR(0,n),    IR(k,k+nb) );

        A1_STAR_MC = A1; 
        B1_VR_STAR = B1;
        Transpose( B1_VR_STAR, B1Trans_STAR_MR, conjugateB );

        // C[MC,MR] += alpha (A1[*,MC])^[T/H] (B1[MR,*])^[T/H]
        //           = alpha (A1^[T/H])[MC,*] (B1^[T/H])[*,MR]
        LocalGemm
        ( orientationOfA, NORMAL, 
          alpha, A1_STAR_MC, B1Trans_STAR_MR, T(1), C );
    }
}

template<typename T>
inline void
SUMMA_TT
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
  T beta,        AbstractDistMatrix<T>& C, GemmAlgorithm alg=GEMM_DEFAULT )
{
    DEBUG_ONLY(CallStackEntry cse("gemm::SUMMA_TT"))
    const Int m = C.Height();
    const Int n = C.Width();
    const Int sumDim = A.Height();
    const double weightTowardsC = 2.;

    switch( alg )
    {
    case GEMM_DEFAULT:
        if( m <= n && weightTowardsC*m <= sumDim )
            SUMMA_TTB( orientationOfA, orientationOfB, alpha, A, B, beta, C );
        else if( n <= m && weightTowardsC*n <= sumDim )
            SUMMA_TTA( orientationOfA, orientationOfB, alpha, A, B, beta, C );
        else
            SUMMA_TTC( orientationOfA, orientationOfB, alpha, A, B, beta, C );
        break;
    case GEMM_SUMMA_A:
        SUMMA_TTA( orientationOfA, orientationOfB, alpha, A, B, beta, C );
        break;
    case GEMM_SUMMA_B:
        SUMMA_TTB( orientationOfA, orientationOfB, alpha, A, B, beta, C );
        break;
    case GEMM_SUMMA_C:
        SUMMA_TTC( orientationOfA, orientationOfB, alpha, A, B, beta, C );
        break;
    default: LogicError("Unsupported Gemm option");
    }
}

} // namespace gemm
} // namespace El
