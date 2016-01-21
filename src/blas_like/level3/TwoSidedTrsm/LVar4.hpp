/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TWOSIDEDTRSM_LVAR4_HPP
#define EL_TWOSIDEDTRSM_LVAR4_HPP

namespace El {
namespace twotrsm {

template<typename F> 
inline void
LVar4( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& L )
{
    DEBUG_ONLY(
      CSE cse("twotrsm::LVar4");
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( L.Height() != L.Width() )
          LogicError("Triangular matrices must be square");
      if( A.Height() != L.Height() )
          LogicError("A and L must be the same size");
    )
    const Int n = A.Height();
    const Int bsize = Blocksize();

    // Temporary products
    Matrix<F> Y21;

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A10 = A( ind1, ind0 );
        auto A11 = A( ind1, ind1 );
        auto A20 = A( ind2, ind0 );
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );

        auto L11 = L( ind1, ind1 );
        auto L21 = L( ind2, ind1 );

        // A10 := inv(L11) A10
        Trsm( LEFT, LOWER, NORMAL, diag, F(1), L11, A10 );

        // A11 := inv(L11) A11 inv(L11)'
        twotrsm::LUnb( diag, A11, L11 );

        // A20 := A20 - L21 A10
        Gemm( NORMAL, NORMAL, F(-1), L21, A10, F(1), A20 );

        // Y21 := L21 A11
        Zeros( Y21, A21.Height(), A21.Width() );
        Hemm( RIGHT, LOWER, F(1), A11, L21, F(0), Y21 );

        // A21 := A21 inv(L11)'
        Trsm( RIGHT, LOWER, ADJOINT, diag, F(1), L11, A21 );

        // A21 := A21 - 1/2 Y21
        Axpy( F(-1)/F(2), Y21, A21 );

        // A22 := A22 - (L21 A21' + A21 L21')
        Her2k( LOWER, NORMAL, F(-1), L21, A21, Base<F>(1), A22 );

        // A21 := A21 - 1/2 Y21
        Axpy( F(-1)/F(2), Y21, A21 );
    }
}

template<typename F> 
inline void
LVar4
( UnitOrNonUnit diag, 
        ElementalMatrix<F>& APre,
  const ElementalMatrix<F>& LPre )
{
    DEBUG_ONLY(
      CSE cse("twotrsm::LVar4");
      if( APre.Height() != APre.Width() )
          LogicError("A must be square");
      if( LPre.Height() != LPre.Width() )
          LogicError("Triangular matrices must be square");
      if( APre.Height() != LPre.Height() )
          LogicError("A and L must be the same size");
    )
    const Int n = APre.Height();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixReadProxy<F,F,MC,MR> LProx( LPre );
    auto& A = AProx.Get();
    auto& L = LProx.GetLocked();

    // Temporary distributions
    DistMatrix<F,STAR,MR  > A10_STAR_MR(g), 
                            A21Adj_STAR_MR(g), L21Adj_STAR_MR(g);
    DistMatrix<F,STAR,MC  > A21Trans_STAR_MC(g);
    DistMatrix<F,STAR,VR  > A10_STAR_VR(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g), L11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g), L21_VC_STAR(g), Y21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> A21_VR_STAR(g), L21_VR_STAR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g);

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A10 = A( ind1, ind0 );
        auto A11 = A( ind1, ind1 );
        auto A20 = A( ind2, ind0 );
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );

        auto L11 = L( ind1, ind1 );
        auto L21 = L( ind2, ind1 );

        // A10 := inv(L11) A10
        L11_STAR_STAR = L11;
        A10_STAR_VR.AlignWith( A20 );
        A10_STAR_VR = A10;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, diag, F(1), L11_STAR_STAR, A10_STAR_VR );

        // A11 := inv(L11) A11 inv(L11)'
        A11_STAR_STAR = A11; 
        TwoSidedTrsm( LOWER, diag, A11_STAR_STAR, L11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A20 := A20 - L21 A10
        L21_MC_STAR.AlignWith( A22 );
        L21_MC_STAR = L21;
        A10_STAR_MR.AlignWith( A20 );
        A10_STAR_MR = A10_STAR_VR;
        LocalGemm( NORMAL, NORMAL, F(-1), L21_MC_STAR, A10_STAR_MR, F(1), A20 );
        A10 = A10_STAR_MR; // delayed write from  A10 := inv(L11) A10

        // Y21 := L21 A11
        L21_VC_STAR.AlignWith( A22 );
        L21_VC_STAR = L21_MC_STAR;
        Y21_VC_STAR.AlignWith( A22 );
        Zeros( Y21_VC_STAR, A21.Height(), nb );
        Hemm
        ( RIGHT, LOWER, 
          F(1), A11_STAR_STAR.Matrix(), L21_VC_STAR.Matrix(), 
          F(0), Y21_VC_STAR.Matrix() );

        // A21 := A21 inv(L11)'
        A21_VC_STAR.AlignWith( A22 );
        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, diag, F(1), L11_STAR_STAR, A21_VC_STAR );

        // A21 := A21 - 1/2 Y21
        Axpy( F(-1)/F(2), Y21_VC_STAR, A21_VC_STAR );

        // A22 := A22 - (L21 A21' + A21 L21')
        A21Trans_STAR_MC.AlignWith( A22 );
        Transpose( A21_VC_STAR, A21Trans_STAR_MC );
        A21_VR_STAR.AlignWith( A22 );
        A21_VR_STAR = A21_VC_STAR;
        L21_VR_STAR.AlignWith( A22 );
        L21_VR_STAR = L21_VC_STAR;
        A21Adj_STAR_MR.AlignWith( A22 );
        L21Adj_STAR_MR.AlignWith( A22 );
        Adjoint( A21_VR_STAR, A21Adj_STAR_MR );
        Adjoint( L21_VR_STAR, L21Adj_STAR_MR );
        LocalTrr2k
        ( LOWER, NORMAL, NORMAL, TRANSPOSE, NORMAL,
          F(-1), L21_MC_STAR,      A21Adj_STAR_MR, 
          F(-1), A21Trans_STAR_MC, L21Adj_STAR_MR,
          F(1),  A22 );

        // A21 := A21 - 1/2 Y21
        Axpy( F(-1)/F(2), Y21_VC_STAR, A21_VC_STAR );
        A21 = A21_VC_STAR;
    }
}

} // namespace twotrsm
} // namespace El

#endif // ifndef EL_TWOSIDEDTRSM_LVAR4_HPP
