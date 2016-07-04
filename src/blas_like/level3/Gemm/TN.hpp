/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace gemm {

// Transpose Normal Gemm that avoids communicating the matrix A
template<typename T> 
void SUMMA_TNA
( Orientation orientA,
  T alpha,
  const AbstractDistMatrix<T>& APre,
  const AbstractDistMatrix<T>& BPre,
        AbstractDistMatrix<T>& CPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( APre, BPre, CPre );
      if( orientA == NORMAL )
          LogicError("A must be (Conjugate)Transposed");
      if( APre.Width() != CPre.Height() || BPre.Width() != CPre.Width() ||
          APre.Height() != BPre.Height() )
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
    DistMatrix<T,MC,STAR> B1_MC_STAR(g);
    DistMatrix<T,MR,STAR> D1_MR_STAR(g);
    DistMatrix<T,MR,MC  > D1_MR_MC(g);

    B1_MC_STAR.AlignWith( A );
    D1_MR_STAR.AlignWith( A );

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        auto B1 = B( ALL, IR(k,k+nb) );
        auto C1 = C( ALL, IR(k,k+nb) );

        // D1[MR,*] := alpha (A1[MC,MR])^T B1[MC,*]
        //           = alpha (A1^T)[MR,MC] B1[MC,*]
        B1_MC_STAR = B1; 
        LocalGemm( orientA, NORMAL, alpha, A, B1_MC_STAR, D1_MR_STAR );

        // C1[MC,MR] += scattered & transposed D1[MR,*] summed over grid cols
        Contract( D1_MR_STAR, D1_MR_MC );
        Axpy( T(1), D1_MR_MC, C1 );
    }
}

// Transpose Normal Gemm that avoids communicating the matrix B
template<typename T> 
void SUMMA_TNB
( Orientation orientA,
  T alpha,
  const AbstractDistMatrix<T>& APre,
  const AbstractDistMatrix<T>& BPre,
        AbstractDistMatrix<T>& CPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( APre, BPre, CPre );
      if( orientA == NORMAL )
          LogicError("A must be (Conjugate)Transposed");
      if( APre.Width() != CPre.Height() || BPre.Width() != CPre.Width() ||
          APre.Height() != BPre.Height() )
          LogicError
          ("Nonconformal matrices:\n",
           DimsString(APre,"A"),"\n",DimsString(BPre,"B"),"\n",
           DimsString(CPre,"C"));
    )
    const Int m = CPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();
    const bool conjugate = ( orientA == ADJOINT );

    DistMatrixReadProxy<T,T,MC,MR> AProx( APre );
    DistMatrixReadProxy<T,T,MC,MR> BProx( BPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> CProx( CPre );
    auto& A = AProx.GetLocked();
    auto& B = BProx.GetLocked();
    auto& C = CProx.Get();

    // Temporary distributions
    DistMatrix<T,MC,STAR> A1_MC_STAR(g);
    DistMatrix<T,MR,STAR> D1Trans_MR_STAR(g);

    A1_MC_STAR.AlignWith( B );
    D1Trans_MR_STAR.AlignWith( B );

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);
        auto A1 = A( ALL,        IR(k,k+nb) );
        auto C1 = C( IR(k,k+nb), ALL        );

        // D1[*,MR] := alpha (A1[MC,*])^[T/H] B[MC,MR]
        //           = alpha (A1^[T/H])[*,MC] B[MC,MR]
        A1_MC_STAR = A1; // A1[MC,*] <- A1[MC,MR]
        LocalGemm( orientA, NORMAL, T(1), B, A1_MC_STAR, D1Trans_MR_STAR );
        TransposeAxpyContract( alpha, D1Trans_MR_STAR, C1, conjugate );
    }
}

// Transpose Normal Gemm that avoids communicating the matrix C
template<typename T> 
void SUMMA_TNC
( Orientation orientA,
  T alpha,
  const AbstractDistMatrix<T>& APre,
  const AbstractDistMatrix<T>& BPre,
        AbstractDistMatrix<T>& CPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( APre, BPre, CPre );
      if( orientA == NORMAL )
          LogicError("A must be (Conjugate)Transposed");
      if( APre.Width() != CPre.Height() || BPre.Width() != CPre.Width() ||
          APre.Height() != BPre.Height() )
          LogicError
          ("Nonconformal matrices:\n",
           DimsString(APre,"A"),"\n",DimsString(BPre,"B"),"\n",
           DimsString(CPre,"C"));
    )
    const Int sumDim = BPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();

    DistMatrixReadProxy<T,T,MC,MR> AProx( APre );
    DistMatrixReadProxy<T,T,MC,MR> BProx( BPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> CProx( CPre );
    auto& A = AProx.GetLocked();
    auto& B = BProx.GetLocked();
    auto& C = CProx.Get();

    // Temporary distributions
    DistMatrix<T,STAR,MC> A1_STAR_MC(g);
    DistMatrix<T,MR,STAR> B1Trans_MR_STAR(g);

    A1_STAR_MC.AlignWith( C );
    B1Trans_MR_STAR.AlignWith( C );

    for( Int k=0; k<sumDim; k+=bsize )
    {
        const Int nb = Min(bsize,sumDim-k);
        auto A1 = A( IR(k,k+nb), ALL );
        auto B1 = B( IR(k,k+nb), ALL );

        // C[MC,MR] += alpha (A1[*,MC])^T B1[*,MR]
        //           = alpha (A1^T)[MC,*] B1[*,MR]
        A1_STAR_MC = A1; 
        Transpose( B1, B1Trans_MR_STAR );
        LocalGemm
        ( orientA, TRANSPOSE, alpha, A1_STAR_MC, B1Trans_MR_STAR, T(1), C );
    }
}

template<typename T>
void SUMMA_TN
( Orientation orientA,
  T alpha,
  const AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& B,
        AbstractDistMatrix<T>& C,
  GemmAlgorithm alg=GEMM_DEFAULT )
{
    DEBUG_CSE
    const Int m = C.Height();
    const Int n = C.Width();
    const Int k = A.Height();
    const double weightTowardsC = 2.;

    switch( alg )
    {
    case GEMM_DEFAULT:
        if( m <= n && weightTowardsC*m <= k )
            SUMMA_TNB( orientA, alpha, A, B, C );
        else if( n <= m && weightTowardsC*n <= k )
            SUMMA_TNA( orientA, alpha, A, B, C );
        else
            SUMMA_TNC( orientA, alpha, A, B, C );
        break;
    case GEMM_SUMMA_A: SUMMA_TNA( orientA, alpha, A, B, C ); break;
    case GEMM_SUMMA_B: SUMMA_TNB( orientA, alpha, A, B, C ); break;
    case GEMM_SUMMA_C: SUMMA_TNC( orientA, alpha, A, B, C ); break;
    default: LogicError("Unsupported Gemm option");
    }
}

} // namespace gemm
} // namespace El
