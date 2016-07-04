/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TWOSIDEDTRMM_LVAR5_HPP
#define EL_TWOSIDEDTRMM_LVAR5_HPP

namespace El {
namespace twotrmm {

// The only reason a field is required is for the existence of 1/2, which is 
// an artifact of the algorithm...
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
    Matrix<F> Y10;

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0, k ),
                         ind1( k, k+nb );

        auto A00 = A( ind0, ind0 );
        auto A10 = A( ind1, ind0 );
        auto A11 = A( ind1, ind1 );

        auto L00 = L( ind0, ind0 );
        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );

        // Y10 := A11 L10
        Y10.Resize( nb, k );
        Zero( Y10 );
        Hemm( LEFT, LOWER, F(1), A11, L10, F(0), Y10 );

        // A10 := A10 L00
        Trmm( RIGHT, LOWER, NORMAL, diag, F(1), L00, A10 );

        // A10 := A10 + 1/2 Y10
        Axpy( F(1)/F(2), Y10, A10 );

        // A00 := A00 + (L10' A10 + A10' L10)
        Her2k( LOWER, ADJOINT, F(1), L10, A10, Base<F>(1), A00 );

        // A10 := A10 + 1/2 Y10
        Axpy( F(1)/F(2), Y10, A10 );

        // A10 := L11' A10
        Trmm( LEFT, LOWER, ADJOINT, diag, F(1), L11, A10 );

        // A11 := L11' A11 L11
        twotrmm::LUnb( diag, A11, L11 );
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
      if( L.Height() != L.Width() )
          LogicError("Triangular matrices must be square");
      if( APre.Height() != L.Height() )
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
    DistMatrix<F,STAR,MC  > A10_STAR_MC(g), L10_STAR_MC(g);
    DistMatrix<F,STAR,MR  > L10_STAR_MR(g);
    DistMatrix<F,STAR,VR  > A10_STAR_VR(g), L10_STAR_VR(g), Y10_STAR_VR(g);
    DistMatrix<F,MR,  STAR> A10Trans_MR_STAR(g);
    DistMatrix<F> Y10(g);

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0, k ),
                         ind1( k, k+nb );

        auto A00 = A( ind0, ind0 );
        auto A10 = A( ind1, ind0 );
        auto A11 = A( ind1, ind1 );

        auto L00 = L( ind0, ind0 );
        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );

        // Y10 := A11 L10
        A11_STAR_STAR = A11;
        L10_STAR_VR.AlignWith( A00 );
        L10_STAR_VR = L10;
        Y10_STAR_VR.AlignWith( A10 );
        Y10_STAR_VR.Resize( nb, k );
        Zero( Y10_STAR_VR );
        Hemm
        ( LEFT, LOWER, 
          F(1), A11_STAR_STAR.Matrix(), L10_STAR_VR.Matrix(),
          F(0), Y10_STAR_VR.Matrix() );
        Y10.AlignWith( A10 );
        Y10 = Y10_STAR_VR;

        // A10 := A10 L00
        Trmm( RIGHT, LOWER, NORMAL, diag, F(1), L00, A10 );

        // A10 := A10 + 1/2 Y10
        Axpy( F(1)/F(2), Y10, A10 );

        // A00 := A00 + (L10' A10 + A10' L10)
        A10Trans_MR_STAR.AlignWith( A00 );
        Transpose( A10, A10Trans_MR_STAR );
        L10_STAR_MR.AlignWith( A00 );
        L10_STAR_MR = L10_STAR_VR;
        A10_STAR_VR.AlignWith( A10 );
        A10_STAR_VR = A10;
        A10_STAR_MC.AlignWith( A00 );
        A10_STAR_MC = A10_STAR_VR;
        L10_STAR_MC.AlignWith( A00 );
        L10_STAR_MC = L10_STAR_VR;
        LocalTrr2k
        ( LOWER, ADJOINT, TRANSPOSE, ADJOINT, NORMAL,
          F(1), L10_STAR_MC, A10Trans_MR_STAR, 
          F(1), A10_STAR_MC, L10_STAR_MR, 
          F(1), A00 );

        // A10 := A10 + 1/2 Y10
        Axpy( F(1)/F(2), Y10_STAR_VR, A10_STAR_VR );

        // A10 := L11' A10
        L11_STAR_STAR = L11;
        LocalTrmm
        ( LEFT, LOWER, ADJOINT, diag, F(1), L11_STAR_STAR, A10_STAR_VR );
        A10 = A10_STAR_VR;

        // A11 := L11' A11 L11
        TwoSidedTrmm( LOWER, diag, A11_STAR_STAR, L11_STAR_STAR );
        A11 = A11_STAR_STAR;
    }
}

} // namespace twotrmm
} // namespace El

#endif // ifndef EL_TWOSIDEDTRMM_LVAR5_HPP
