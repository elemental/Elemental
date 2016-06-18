/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_INVERSE_LUPARTIALPIV_HPP
#define EL_INVERSE_LUPARTIALPIV_HPP

namespace El {
namespace inverse {

// Start by forming the partially pivoted LU decomposition of A,
//     P A = L U,
// then inverting the system of equations,
//     inv(A) inv(P) = inv(U) inv(L),
// then,
//     inv(A) = inv(U) inv(L) P.

template<typename F> 
void AfterLUPartialPiv( Matrix<F>& A, const Permutation& P )
{
    DEBUG_CSE
    if( A.Height() != A.Width() )
        LogicError("Cannot invert non-square matrices");

    TriangularInverse( UPPER, NON_UNIT, A );

    const Int n = A.Height();

    // Solve inv(A) L = inv(U) for inv(A)
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);
        const IR ind1( k, k+nb ), ind2( k+nb, END );

        auto A1 = A( ALL, ind1 );
        auto A2 = A( ALL, ind2 );

        auto A11 = A( ind1, ind1 );
        auto A21 = A( ind2, ind1 );

        // Copy out L1
        auto L11( A11 );
        auto L21( A21 );

        // Zero the strictly lower triangular portion of A1
        MakeTrapezoidal( UPPER, A11 );
        Zero( A21 );

        // Perform the lazy update of A1
        Gemm( NORMAL, NORMAL, F(-1), A2, L21, F(1), A1 );

        // Solve against this diagonal block of L11
        Trsm( RIGHT, LOWER, NORMAL, UNIT, F(1), L11, A1 );
    }

    // inv(A) := inv(A) P
    P.InversePermuteCols( A );
}

template<typename F> 
void LUPartialPiv( Matrix<F>& A )
{
    DEBUG_CSE
    if( A.Height() != A.Width() )
        LogicError("Cannot invert non-square matrices");
    Permutation P;
    El::LU( A, P );
    inverse::AfterLUPartialPiv( A, P );
}

template<typename F> 
void AfterLUPartialPiv
(       ElementalMatrix<F>& APre,
  const DistPermutation& P )
{
    DEBUG_CSE

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    if( A.Height() != A.Width() )
        LogicError("Cannot invert non-square matrices");

    TriangularInverse( UPPER, NON_UNIT, A );

    const Grid& g = A.Grid();
    DistMatrix<F,VC,  STAR> A1_VC_STAR(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,VR,  STAR> L21_VR_STAR(g);
    DistMatrix<F,STAR,MR  > L21Trans_STAR_MR(g);
    DistMatrix<F,MC,  STAR> Z1(g);

    const Int n = A.Height();

    // Solve inv(A) L = inv(U) for inv(A)
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);
        const IR ind1( k, k+nb ), ind2( k+nb, END );

        auto A1 = A( ALL, ind1 );
        auto A2 = A( ALL, ind2 );

        auto A11 = A( ind1, ind1 );
        auto A21 = A( ind2, ind1 );

        // Copy out L1
        L11_STAR_STAR = A11;
        L21_VR_STAR.AlignWith( A2 );
        L21_VR_STAR = A21;
        L21Trans_STAR_MR.AlignWith( A2 );
        Transpose( L21_VR_STAR, L21Trans_STAR_MR );

        // Zero the strictly lower triangular portion of A1
        MakeTrapezoidal( UPPER, A11 );
        Zero( A21 );

        // Perform the lazy update of A1
        Z1.AlignWith( A1 );
        Zeros( Z1, n, nb );
        LocalGemm( NORMAL, TRANSPOSE, F(-1), A2, L21Trans_STAR_MR, F(0), Z1 );
        AxpyContract( F(1), Z1, A1 );

        // Solve against this diagonal block of L11
        A1_VC_STAR = A1;
        LocalTrsm
        ( RIGHT, LOWER, NORMAL, UNIT, F(1), L11_STAR_STAR, A1_VC_STAR );
        A1 = A1_VC_STAR;
    }

    // inv(A) := inv(A) P
    P.InversePermuteCols( A );
}

template<typename F> 
void LUPartialPiv( ElementalMatrix<F>& A )
{
    DEBUG_CSE
    if( A.Height() != A.Width() )
        LogicError("Cannot invert non-square matrices");
    const Grid& g = A.Grid();
    DistPermutation P(g);
    El::LU( A, P );
    inverse::AfterLUPartialPiv( A, P );
}

} // namespace inverse
} // namespace El

#endif // ifndef EL_INVERSE_LUPARTIALPIV_HPP
