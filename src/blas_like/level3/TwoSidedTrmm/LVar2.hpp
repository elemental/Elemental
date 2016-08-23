/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TWOSIDEDTRMM_LVAR2_HPP
#define EL_TWOSIDEDTRMM_LVAR2_HPP

namespace El {
namespace twotrmm {

// The only reason a field is required is for the existence of 1/2, which is 
// an artifact of the algorithm...
template<typename F> 
void LVar2( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& L )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError( "A must be square." );
      if( L.Height() != L.Width() )
          LogicError( "Triangular matrices must be square." );
      if( A.Height() != L.Height() )
          LogicError( "A and L must be the same size." );
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
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );

        auto L11 = L( ind1, ind1 );
        auto L21 = L( ind2, ind1 );

        // A10 := L11' A10
        Trmm( LEFT, LOWER, ADJOINT, diag, F(1), L11, A10 );

        // A10 := A10 + L21' A20
        Gemm( ADJOINT, NORMAL, F(1), L21, A20, F(1), A10 );

        // Y21 := A22 L21
        Y21.Resize( A21.Height(), nb );
        Zero( Y21 );
        Hemm( LEFT, LOWER, F(1), A22, L21, F(0), Y21 );

        // A21 := A21 L11
        Trmm( RIGHT, LOWER, NORMAL, diag, F(1), L11, A21 );

        // A21 := A21 + 1/2 Y21
        Axpy( F(1)/F(2), Y21, A21 );

        // A11 := L11' A11 L11
        twotrmm::LUnb( diag, A11, L11 );

        // A11 := A11 + (A21' L21 + L21' A21)
        Her2k( LOWER, ADJOINT, Base<F>(1), A21, L21, Base<F>(1), A11 );

        // A21 := A21 + 1/2 Y21
        Axpy( F(1)/F(2), Y21, A21 );
    }
}

template<typename F> 
void LVar2
( UnitOrNonUnit diag, 
        AbstractDistMatrix<F>& APre,
  const AbstractDistMatrix<F>& LPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( APre.Height() != APre.Width() )
          LogicError( "A must be square." );
      if( LPre.Height() != LPre.Width() )
          LogicError( "Triangular matrices must be square." );
      if( APre.Height() != LPre.Height() )
          LogicError( "A and L must be the same size." );
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
    DistMatrix<F,STAR,MR  > L21Adj_STAR_MR(g), X10_STAR_MR(g);
    DistMatrix<F,STAR,VR  > A10_STAR_VR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g), Z21_MC_STAR(g);
    DistMatrix<F,MR,  STAR> Z21_MR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g), L21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> L21_VR_STAR(g);
    DistMatrix<F,MR,  MC  > Z21_MR_MC(g);
    DistMatrix<F> Y21(g);

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A10 = A( ind1, ind0 );
        auto A11 = A( ind1, ind1 );
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );

        auto L11 = L( ind1, ind1 );
        auto L21 = L( ind2, ind1 );

        // A10 := L11' A10
        L11_STAR_STAR = L11;
        A10_STAR_VR = A10;
        LocalTrmm
        ( LEFT, LOWER, ADJOINT, diag, F(1), L11_STAR_STAR, A10_STAR_VR );
        A10 = A10_STAR_VR;

        // A10 := A10 + L21' A20
        L21_MC_STAR.AlignWith( A20 );
        L21_MC_STAR = L21;
        X10_STAR_MR.AlignWith( A10 );
        LocalGemm( ADJOINT, NORMAL, F(1), L21_MC_STAR, A20, X10_STAR_MR );
        AxpyContract( F(1), X10_STAR_MR, A10 );

        // Y21 := A22 L21
        L21_VC_STAR.AlignWith( A22 );
        L21_VR_STAR.AlignWith( A22 );
        L21_VC_STAR = L21_MC_STAR;
        L21_VR_STAR = L21_VC_STAR;
        L21Adj_STAR_MR.AlignWith( A22 );
        L21_VR_STAR.AdjointPartialColAllGather( L21Adj_STAR_MR );
        Z21_MC_STAR.AlignWith( A22 );
        Z21_MR_STAR.AlignWith( A22 );
        Z21_MC_STAR.Resize( A21.Height(), nb );
        Z21_MR_STAR.Resize( A21.Height(), nb );
        Zero( Z21_MC_STAR );
        Zero( Z21_MR_STAR );
        symm::LocalAccumulateLL
        ( ADJOINT, 
          F(1), A22, L21_MC_STAR, L21Adj_STAR_MR, Z21_MC_STAR, Z21_MR_STAR );
        Contract( Z21_MR_STAR, Z21_MR_MC );
        Y21.AlignWith( A21 );
        Y21 = Z21_MR_MC;
        AxpyContract( F(1), Z21_MC_STAR, Y21 );

        // A21 := A21 L11
        A21_VC_STAR.AlignWith( A22 );
        A21_VC_STAR = A21;
        LocalTrmm
        ( RIGHT, LOWER, NORMAL, diag, F(1), L11_STAR_STAR, A21_VC_STAR );
        A21 = A21_VC_STAR;

        // A21 := A21 + 1/2 Y21
        Axpy( F(1)/F(2), Y21, A21 );

        // A11 := L11' A11 L11
        A11_STAR_STAR = A11;
        TwoSidedTrmm( LOWER, diag, A11_STAR_STAR, L11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A11 := A11 + (A21' L21 + L21' A21)
        A21_VC_STAR = A21;
        X11_STAR_STAR.Resize( nb, nb );
        Zero( X11_STAR_STAR );
        Her2k
        ( LOWER, ADJOINT,
          F(1), A21_VC_STAR.Matrix(), L21_VC_STAR.Matrix(),
          F(0), X11_STAR_STAR.Matrix() );
        AxpyContract( F(1), X11_STAR_STAR, A11 );

        // A21 := A21 + 1/2 Y21
        Axpy( F(1)/F(2), Y21, A21 );
    }
}

} // namespace twotrmm
} // namespace El

#endif // ifndef EL_TWOSIDEDTRMM_LVAR2_HPP
