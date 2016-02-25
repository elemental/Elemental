/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TWOSIDEDTRSM_LVAR3_HPP
#define EL_TWOSIDEDTRSM_LVAR3_HPP

namespace El {
namespace twotrsm {

template<typename F>
inline void
LVar3( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& L )
{
    DEBUG_ONLY(
      CSE cse("twotrsm::LVar3");
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( L.Height() != L.Width() )
          LogicError("Triangular matrices must be square");
      if( A.Height() != L.Height() )
          LogicError("A and L must be the same size");
    )
    const Int n = A.Height();
    const Int bsize = Blocksize();

    // We will use an entire extra matrix as temporary storage. If this is not
    // acceptable, use LVar4 instead.
    Matrix<F> Y;
    Zeros( Y, n, n );

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

        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );
        auto L21 = L( ind2, ind1 );

        auto Y10 = Y( ind1, ind0 );
        auto Y20 = Y( ind2, ind0 );
        auto Y21 = Y( ind2, ind1 );

        // A10 := A10 - 1/2 Y10
        Axpy( F(-1)/F(2), Y10, A10 );

        // A11 := A11 - (A10 L10' + L10 A10')
        Her2k( LOWER, NORMAL, F(-1), A10, L10, F(1), A11 );

        // A11 := inv(L11) A11 inv(L11)'
        twotrsm::LUnb( diag, A11, L11 );

        // A21 := A21 - A20 L10'
        Gemm( NORMAL, ADJOINT, F(-1), A20, L10, F(1), A21 );

        // A21 := A21 inv(L11)'
        Trsm( RIGHT, LOWER, ADJOINT, diag, F(1), L11, A21 );

        // A10 := A10 - 1/2 Y10
        Axpy( F(-1)/F(2), Y10, A10 );

        // A10 := inv(L11) A10
        Trsm( LEFT, LOWER, NORMAL, diag, F(1), L11, A10 );

        // Y20 := Y20 + L21 A10
        Gemm( NORMAL, NORMAL, F(1), L21, A10, F(1), Y20 );

        // Y21 := L21 A11
        Hemm( RIGHT, LOWER, F(1), A11, L21, F(0), Y21 );

        // Y21 := Y21 + L20 A10'
        Gemm( NORMAL, ADJOINT, F(1), L20, A10, F(1), Y21 );
    }
}

template<typename F>
inline void
LVar3
( UnitOrNonUnit diag, 
        ElementalMatrix<F>& APre,
  const ElementalMatrix<F>& LPre )
{
    DEBUG_ONLY(
      CSE cse("twotrsm::LVar3");
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
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g), L11_STAR_STAR(g), 
                            X11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,STAR,VR  > A10_STAR_VR(g), L10_STAR_VR(g);
    DistMatrix<F,STAR,MR  > A10_STAR_MR(g), L10_STAR_MR(g), A11_STAR_MR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g), X21_MC_STAR(g), Z21_MC_STAR(g);

    // We will use an entire extra matrix as temporary storage. If this is not
    // acceptable, use LVar4 instead.
    DistMatrix<F> Y(g);
    Y.AlignWith( A );
    Zeros( Y, n, n );

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

        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );
        auto L21 = L( ind2, ind1 );

        auto Y10 = Y( ind1, ind0 );
        auto Y20 = Y( ind2, ind0 );
        auto Y21 = Y( ind2, ind1 );

        // A10 := A10 - 1/2 Y10
        Axpy( F(-1)/F(2), Y10, A10 );

        // A11 := A11 - (A10 L10' + L10 A10')
        A10_STAR_VR.AlignWith( A10 );
        L10_STAR_VR.AlignWith( A10 );
        A10_STAR_VR = A10;
        L10_STAR_VR = L10;
        Zeros( X11_STAR_STAR, nb, nb );
        Her2k
        ( LOWER, NORMAL, 
          F(1), A10_STAR_VR.Matrix(), L10_STAR_VR.Matrix(),
          F(0), X11_STAR_STAR.Matrix() );
        MakeTrapezoidal( LOWER, X11_STAR_STAR );
        AxpyContract( F(-1), X11_STAR_STAR, A11 );

        // A11 := inv(L11) A11 inv(L11)'
        A11_STAR_STAR = A11;
        L11_STAR_STAR = L11;
        TwoSidedTrsm( LOWER, diag, A11_STAR_STAR, L11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A21 := A21 - A20 L10'
        L10_STAR_MR.AlignWith( A10 );
        L10_STAR_MR = L10_STAR_VR;
        X21_MC_STAR.AlignWith( A20 );
        LocalGemm( NORMAL, ADJOINT, F(1), A20, L10_STAR_MR, X21_MC_STAR );
        AxpyContract( F(-1), X21_MC_STAR, A21 );

        // A21 := A21 inv(L11)'
        A21_VC_STAR.AlignWith( A21 );
        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, diag, F(1), L11_STAR_STAR, A21_VC_STAR );
        A21 = A21_VC_STAR;

        // A10 := A10 - 1/2 Y10
        Axpy( F(-1)/F(2), Y10, A10 );

        // A10 := inv(L11) A10
        A10_STAR_VR = A10;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, diag, F(1), L11_STAR_STAR, A10_STAR_VR );

        // Y20 := Y20 + L21 A10
        A10_STAR_MR.AlignWith( A10 );
        A10_STAR_MR = A10_STAR_VR;
        A10 = A10_STAR_MR;
        L21_MC_STAR.AlignWith( Y21 );
        L21_MC_STAR = L21;
        LocalGemm( NORMAL, NORMAL, F(1), L21_MC_STAR, A10_STAR_MR, F(1), Y20 );

        // Y21 := L21 A11
        MakeHermitian( LOWER, A11_STAR_STAR );
        A11_STAR_MR.AlignWith( Y21 );
        A11_STAR_MR = A11_STAR_STAR;
        LocalGemm( NORMAL, NORMAL, F(1), L21_MC_STAR, A11_STAR_MR, F(0), Y21 );

        // Y21 := Y21 + L20 A10'
        Z21_MC_STAR.AlignWith( L20 );
        LocalGemm( NORMAL, ADJOINT, F(1), L20, A10_STAR_MR, Z21_MC_STAR );
        AxpyContract( F(1), Z21_MC_STAR, Y21 );
    }
}

} // namespace twotrsm
} // namespace El

#endif // ifndef EL_TWOSIDEDTRSM_LVAR3_HPP
