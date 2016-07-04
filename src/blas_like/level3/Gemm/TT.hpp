/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace gemm {

// Transpose Transpose Gemm that avoids communicating the matrix A
template<typename T>
void SUMMA_TTA
( Orientation orientA,
  Orientation orientB,
  T alpha,
  const AbstractDistMatrix<T>& APre,
  const AbstractDistMatrix<T>& BPre,
        AbstractDistMatrix<T>& CPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( APre, BPre, CPre );
      if( orientA == NORMAL || orientB == NORMAL )
          LogicError("A and B must be (Conjugate)Transposed");
      if( APre.Width() != CPre.Height() || BPre.Height() != CPre.Width() ||
          APre.Height() != BPre.Width() )
          LogicError
          ("Nonconformal matrices:\n",
           DimsString(APre,"A"),"\n",DimsString(BPre,"B"),"\n",
           DimsString(CPre,"C"));
    )
    const Int n = CPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();

    DistMatrixReadProxy<T,T,MC,MR> AProx( APre );
    DistMatrixReadProxy<T,T,MC,MR> BProx( BPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> CProx( CPre );
    auto& A = AProx.GetLocked();
    auto& B = BProx.GetLocked();
    auto& C = CProx.Get();

    // Temporary distributions
    DistMatrix<T,STAR,MC  > B1_STAR_MC(g);
    DistMatrix<T,MR,  MC  > D1_MR_MC(g);
    DistMatrix<T,MR,  STAR> D1_MR_STAR(g);

    B1_STAR_MC.AlignWith( A ); 
    D1_MR_STAR.AlignWith( A );  

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        auto B1 = B( IR(k,k+nb), ALL        );
        auto C1 = C( ALL,        IR(k,k+nb) );

        B1_STAR_MC = B1; // B1[*,MC] <- B1[MC,MR]

        // D1[MR,*] := alpha (A[MC,MR])^T (B1[*,MC])^T
        //           = alpha (A^T)[MR,MC] (B1^T)[MC,*]
        LocalGemm( orientA, orientB, alpha, A, B1_STAR_MC, D1_MR_STAR );

        // C1[MC,MR] += scattered D1[MR,*] summed over grid cols
        Contract( D1_MR_STAR, D1_MR_MC );
        Axpy( T(1), D1_MR_MC, C1 );
    }
}

// Transpose Transpose Gemm that avoids communicating the matrix B
template<typename T>
void SUMMA_TTB
( Orientation orientA,
  Orientation orientB,
  T alpha,
  const AbstractDistMatrix<T>& APre,
  const AbstractDistMatrix<T>& BPre,
        AbstractDistMatrix<T>& CPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( APre, BPre, CPre );
      if( orientA == NORMAL || orientB == NORMAL )
          LogicError("A and B must be (Conjugate)Transposed");
      if( APre.Width() != CPre.Height() || BPre.Height() != CPre.Width() ||
          APre.Height() != BPre.Width() )
          LogicError
          ("Nonconformal matrices:\n",
           DimsString(APre,"A"),"\n",DimsString(BPre,"B"),"\n",
           DimsString(CPre,"C"));
    )
    const Int m = CPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();
    const bool conjugateA = ( orientA == ADJOINT ); 

    DistMatrixReadProxy<T,T,MC,MR> AProx( APre );
    DistMatrixReadProxy<T,T,MC,MR> BProx( BPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> CProx( CPre );
    auto& A = AProx.GetLocked();
    auto& B = BProx.GetLocked();
    auto& C = CProx.Get();

    // Temporary distributions
    DistMatrix<T,VR,  STAR> A1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > A1Trans_STAR_MR(g);
    DistMatrix<T,STAR,MC  > D1_STAR_MC(g);
    DistMatrix<T,MR,  MC  > D1_MR_MC(g);

    A1_VR_STAR.AlignWith( B );
    A1Trans_STAR_MR.AlignWith( B );
    D1_STAR_MC.AlignWith( B );

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);
        auto A1 = A( ALL,        IR(k,k+nb) );
        auto C1 = C( IR(k,k+nb), ALL        );

        // D1[*,MC] := alpha (A1[MR,*])^[T/H] (B[MC,MR])^[T/H]
        //           = alpha (A1^[T/H])[*,MR] (B^[T/H])[MR,MC]
        A1_VR_STAR = A1;
        Transpose( A1_VR_STAR, A1Trans_STAR_MR, conjugateA );
        LocalGemm( NORMAL, orientB, alpha, A1Trans_STAR_MR, B, D1_STAR_MC );

        // C1[MC,MR] += scattered & transposed D1[*,MC] summed over grid rows
        Contract( D1_STAR_MC, D1_MR_MC );
        Axpy( T(1), D1_MR_MC, C1 );
    }
}

// Transpose Transpose Gemm that avoids communicating the matrix C
template<typename T>
void SUMMA_TTC
( Orientation orientA,
  Orientation orientB,
  T alpha,
  const AbstractDistMatrix<T>& APre,
  const AbstractDistMatrix<T>& BPre,
        AbstractDistMatrix<T>& CPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( APre, BPre, CPre );
      if( orientA == NORMAL || orientB == NORMAL )
          LogicError("A and B must be (Conjugate)Transposed");
      if( APre.Width() != CPre.Height() || BPre.Height() != CPre.Width() ||
          APre.Height() != BPre.Width() )
          LogicError
          ("Nonconformal matrices:\n",
           DimsString(APre,"A"),"\n",DimsString(BPre,"B"),"\n",
           DimsString(CPre,"C"));
    )
    const Int sumDim = APre.Height();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();
    const bool conjugateB = ( orientB == ADJOINT );

    DistMatrixReadProxy<T,T,MC,MR> AProx( APre );
    DistMatrixReadProxy<T,T,MC,MR> BProx( BPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> CProx( CPre );
    auto& A = AProx.GetLocked();
    auto& B = BProx.GetLocked();
    auto& C = CProx.Get();

    // Temporary distributions
    DistMatrix<T,STAR,MC  > A1_STAR_MC(g);
    DistMatrix<T,VR,  STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > B1Trans_STAR_MR(g);

    A1_STAR_MC.AlignWith( C );
    B1_VR_STAR.AlignWith( C );
    B1Trans_STAR_MR.AlignWith( C );
    
    for( Int k=0; k<sumDim; k+=bsize )
    {
        const Int nb = Min(bsize,sumDim-k);
        auto A1 = A( IR(k,k+nb), ALL        );
        auto B1 = B( ALL,        IR(k,k+nb) );

        A1_STAR_MC = A1; 
        B1_VR_STAR = B1;
        Transpose( B1_VR_STAR, B1Trans_STAR_MR, conjugateB );

        // C[MC,MR] += alpha (A1[*,MC])^[T/H] (B1[MR,*])^[T/H]
        //           = alpha (A1^[T/H])[MC,*] (B1^[T/H])[*,MR]
        LocalGemm
        ( orientA, NORMAL, alpha, A1_STAR_MC, B1Trans_STAR_MR, T(1), C );
    }
}

template<typename T>
void SUMMA_TT
( Orientation orientA,
  Orientation orientB,
  T alpha,
  const AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& B,
        AbstractDistMatrix<T>& C,
  GemmAlgorithm alg=GEMM_DEFAULT )
{
    DEBUG_CSE
    const Int m = C.Height();
    const Int n = C.Width();
    const Int sumDim = A.Height();
    const double weightTowardsC = 2.;

    switch( alg )
    {
    case GEMM_DEFAULT:
        if( m <= n && weightTowardsC*m <= sumDim )
            SUMMA_TTB( orientA, orientB, alpha, A, B, C );
        else if( n <= m && weightTowardsC*n <= sumDim )
            SUMMA_TTA( orientA, orientB, alpha, A, B, C );
        else
            SUMMA_TTC( orientA, orientB, alpha, A, B, C );
        break;
    case GEMM_SUMMA_A:
        SUMMA_TTA( orientA, orientB, alpha, A, B, C );
        break;
    case GEMM_SUMMA_B:
        SUMMA_TTB( orientA, orientB, alpha, A, B, C );
        break;
    case GEMM_SUMMA_C:
        SUMMA_TTC( orientA, orientB, alpha, A, B, C );
        break;
    default: LogicError("Unsupported Gemm option");
    }
}

} // namespace gemm
} // namespace El
