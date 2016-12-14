/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TWOSIDEDTRSM_LVAR1_HPP
#define EL_TWOSIDEDTRSM_LVAR1_HPP

namespace El {
namespace twotrsm {

template<typename F> 
void LVar1( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& L )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
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

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto A00 = A( ind0, ind0 );
        auto A10 = A( ind1, ind0 );
        auto A11 = A( ind1, ind1 );

        auto L00 = L( ind0, ind0 );
        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );

        // Y10 := L10 A00
        L10.Resize( L10.Height(), A00.Width() );
        Zero( L10 );
        Hemm( RIGHT, LOWER, F(1), A00, L10, F(0), Y10 );

        // A10 := A10 inv(L00)'
        Trsm( RIGHT, LOWER, ADJOINT, diag, F(1), L00, A10 );

        // A10 := A10 - 1/2 Y10
        Axpy( F(-1)/F(2), Y10, A10 );

        // A11 := A11 - (A10 L10' + L10 A10')
        Her2k( LOWER, NORMAL, F(-1), A10, L10, F(1), A11 );

        // A11 := inv(L11) A11 inv(L11)'
        twotrsm::LUnb( diag, A11, L11 );

        // A10 := A10 - 1/2 Y10
        Axpy( F(-1)/F(2), Y10, A10 );

        // A10 := inv(L11) A10
        Trsm( LEFT, LOWER, NORMAL, diag, F(1), L11, A10 );
    }
}

template<typename F> 
void LVar1
( UnitOrNonUnit diag, 
        AbstractDistMatrix<F>& APre,
  const AbstractDistMatrix<F>& LPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
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
    DistMatrix<F,STAR,VR  > A10_STAR_VR(g), L10_STAR_VR(g);
    DistMatrix<F,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<F,MR,  STAR> L10Adj_MR_STAR(g), Z10Adj_MR_STAR(g);
    DistMatrix<F,MC,  STAR> Z10Adj_MC_STAR(g);
    DistMatrix<F,VC,  STAR> L10Adj_VC_STAR(g);
    DistMatrix<F,MR,  MC  > Z10Adj_MR_MC(g);
    DistMatrix<F> Y10(g), Z10Adj(g);

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto A00 = A( ind0, ind0 );
        auto A10 = A( ind1, ind0 );
        auto A11 = A( ind1, ind1 );

        auto L00 = L( ind0, ind0 );
        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );

        // Y10 := L10 A00
        L10Adj_MR_STAR.AlignWith( A00 );
        L10.AdjointColAllGather( L10Adj_MR_STAR );
        L10Adj_VC_STAR.AlignWith( A00 );
        L10Adj_VC_STAR = L10Adj_MR_STAR;
        L10_STAR_MC.AlignWith( A00 );
        L10Adj_VC_STAR.AdjointPartialColAllGather( L10_STAR_MC );
        Z10Adj_MC_STAR.AlignWith( A00 );
        Z10Adj_MR_STAR.AlignWith( A00 );
        Z10Adj_MC_STAR.Resize( k, nb );
        Z10Adj_MR_STAR.Resize( k, nb );
        Zero( Z10Adj_MC_STAR );
        Zero( Z10Adj_MR_STAR );
        symm::LocalAccumulateRL
        ( ADJOINT,
          F(1), A00, L10_STAR_MC, L10Adj_MR_STAR, 
          Z10Adj_MC_STAR, Z10Adj_MR_STAR );
        Z10Adj.AlignWith( A10 );
        Contract( Z10Adj_MC_STAR, Z10Adj );
        Z10Adj_MR_MC.AlignWith( A10 );
        Z10Adj_MR_MC = Z10Adj;
        AxpyContract( F(1), Z10Adj_MR_STAR, Z10Adj_MR_MC );
        Y10.AlignWith( A10 );
        Adjoint( Z10Adj_MR_MC, Y10 );

        // A10 := A10 inv(L00)'
        // This is the bottleneck because A10 only has blocksize rows
        Trsm( RIGHT, LOWER, ADJOINT, diag, F(1), L00, A10 );

        // A10 := A10 - 1/2 Y10
        Axpy( F(-1)/F(2), Y10, A10 );

        // A11 := A11 - (A10 L10' + L10 A10')
        A10_STAR_VR.AlignWith( A10 );
        A10_STAR_VR = A10;
        L10_STAR_VR.AlignWith( A00 );
        L10_STAR_VR = L10;
        X11_STAR_STAR.Resize( nb, nb );
        Zero( X11_STAR_STAR );
        Her2k
        ( LOWER, NORMAL,
          F(-1), A10_STAR_VR.Matrix(), L10_STAR_VR.Matrix(), 
          F(0), X11_STAR_STAR.Matrix() );
        AxpyContract( F(1), X11_STAR_STAR, A11 );

        // A11 := inv(L11) A11 inv(L11)'
        A11_STAR_STAR = A11;
        L11_STAR_STAR = L11;
        TwoSidedTrsm( LOWER, diag, A11_STAR_STAR, L11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A10 := A10 - 1/2 Y10
        Axpy( F(-1)/F(2), Y10, A10 );

        // A10 := inv(L11) A10
        A10_STAR_VR = A10;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, diag, F(1), L11_STAR_STAR, A10_STAR_VR );
        A10 = A10_STAR_VR;
    }
}

} // namespace twotrsm
} // namespace El

#endif // ifndef EL_TWOSIDEDTRSM_LVAR1_HPP
