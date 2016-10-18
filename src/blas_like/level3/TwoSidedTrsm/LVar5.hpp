/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TWOSIDEDTRSM_LVAR5_HPP
#define EL_TWOSIDEDTRSM_LVAR5_HPP

namespace El {
namespace twotrsm {

template<typename F> 
void LVar5( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& L )
{
    DEBUG_CSE
    DEBUG_ONLY(
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

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A11 = A( ind1, ind1 );
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );

        auto L11 = L( ind1, ind1 );
        auto L21 = L( ind2, ind1 );

        // A11 := inv(L11) A11 inv(L11)'
        twotrsm::LUnb( diag, A11, L11 );

        // Y21 := L21 A11
        Y21.Resize( A21.Height(), A21.Width() );
        Zero( Y21 );
        Hemm( RIGHT, LOWER, F(1), A11, L21, F(0), Y21 );

        // A21 := A21 inv(L11)'
        Trsm( RIGHT, LOWER, ADJOINT, diag, F(1), L11, A21 );

        // A21 := A21 - 1/2 Y21
        Axpy( F(-1)/F(2), Y21, A21 );

        // A22 := A22 - (L21 A21' + A21 L21')
        Her2k( LOWER, NORMAL, F(-1), L21, A21, Base<F>(1), A22 );

        // A21 := A21 - 1/2 Y21
        Axpy( F(-1)/F(2), Y21, A21 );

        // A21 := inv(L22) A21
        Trsm( LEFT, LOWER, NORMAL, diag, F(1), L22, A21 );
    }
}

template<typename F> 
void LVar5
( UnitOrNonUnit diag, 
        AbstractDistMatrix<F>& APre,
  const AbstractDistMatrix<F>& LPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g), L11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> A21_MC_STAR(g), L21_MC_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g), L21_VC_STAR(g), Y21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> A21_VR_STAR(g), L21_VR_STAR(g);
    DistMatrix<F,STAR,MR  > A21Adj_STAR_MR(g), L21Adj_STAR_MR(g);
    DistMatrix<F> Y21(g);

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A11 = A( ind1, ind1 );
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );

        auto L11 = L( ind1, ind1 );
        auto L21 = L( ind2, ind1 );

        // A11 := inv(L11) A11 inv(L11)'
        L11_STAR_STAR = L11;
        A11_STAR_STAR = A11;
        TwoSidedTrsm( LOWER, diag, A11_STAR_STAR, L11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // Y21 := L21 A11
        L21_VC_STAR.AlignWith( A22 );
        L21_VC_STAR = L21;
        Y21_VC_STAR.AlignWith( A22 );
        Y21_VC_STAR.Resize( A21.Height(), A21.Width() );
        Zero( Y21_VC_STAR );
        Hemm
        ( RIGHT, LOWER, 
          F(1), A11_STAR_STAR.Matrix(), L21_VC_STAR.Matrix(), 
          F(0), Y21_VC_STAR.Matrix() );
        Y21.AlignWith( A21 );
        Y21 = Y21_VC_STAR;

        // A21 := A21 inv(L11)'
        A21_VC_STAR.AlignWith( A22 );
        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, diag, F(1), L11_STAR_STAR, A21_VC_STAR );
        A21 = A21_VC_STAR;

        // A21 := A21 - 1/2 Y21
        Axpy( F(-1)/F(2), Y21, A21 );

        // A22 := A22 - (L21 A21' + A21 L21')
        A21_MC_STAR.AlignWith( A22 );
        L21_MC_STAR.AlignWith( A22 );
        A21_MC_STAR = A21;
        L21_MC_STAR = L21;
        A21_VC_STAR = A21_MC_STAR;
        A21_VR_STAR.AlignWith( A22 );
        L21_VR_STAR.AlignWith( A22 );
        A21_VR_STAR = A21_VC_STAR;
        L21_VR_STAR = L21_VC_STAR;
        A21Adj_STAR_MR.AlignWith( A22 );
        L21Adj_STAR_MR.AlignWith( A22 );
        A21_VR_STAR.AdjointPartialColAllGather( A21Adj_STAR_MR );
        L21_VR_STAR.AdjointPartialColAllGather( L21Adj_STAR_MR );
        LocalTrr2k
        ( LOWER, NORMAL, NORMAL, NORMAL, NORMAL,
          F(-1), L21_MC_STAR, A21Adj_STAR_MR,
          F(-1), A21_MC_STAR, L21Adj_STAR_MR, F(1), A22 );

        // A21 := A21 - 1/2 Y21
        Axpy( F(-1)/F(2), Y21, A21 );

        // A21 := inv(L22) A21
        //
        // This is the bottleneck because A21 only has blocksize columns
        Trsm( LEFT, LOWER, NORMAL, diag, F(1), L22, A21 );
    }
}

} // namespace twotrsm
} // namespace El

#endif // ifndef EL_TWOSIDEDTRSM_LVAR5_HPP
