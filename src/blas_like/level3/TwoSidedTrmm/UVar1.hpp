/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TWOSIDEDTRMM_UVAR1_HPP
#define EL_TWOSIDEDTRMM_UVAR1_HPP

namespace El {
namespace twotrmm {

// The only reason a field is required is for the existence of 1/2, which is 
// an artifact of the algorithm...
template<typename F> 
void UVar1( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& U )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( U.Height() != U.Width() )
          LogicError("Triangular matrices must be square");
      if( A.Height() != U.Height() )
          LogicError("A and U must be the same size");
    )
    const Int n = A.Height();
    const Int bsize = Blocksize();

    // Temporary products
    Matrix<F> Y12;

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );
        auto A22 = A( ind2, ind2 );

        auto U11 = U( ind1, ind1 );
        auto U12 = U( ind1, ind2 );
        auto U22 = U( ind2, ind2 );

        // Y12 := U12 A22
        Y12.Resize( nb, A12.Width() );
        Zero( Y12 );
        Hemm( RIGHT, UPPER, F(1), A22, U12, F(0), Y12 );

        // A12 := U11 A12
        Trmm( LEFT, UPPER, NORMAL, diag, F(1), U11, A12 );

        // A12 := A12 + 1/2 Y12
        Axpy( F(1)/F(2), Y12, A12 );

        // A11 := U11 A11 U11'
        twotrmm::UUnb( diag, A11, U11 );

        // A11 := A11 + (U12 A12' + A12 U12')
        Her2k( UPPER, NORMAL, Base<F>(1), U12, A12, Base<F>(1), A11 );

        // A12 := A12 + 1/2 Y12
        Axpy( F(1)/F(2), Y12, A12 );

        // A12 := A12 U22'
        Trmm( RIGHT, UPPER, ADJOINT, diag, F(1), U22, A12 );
    }
}

template<typename F> 
void UVar1
( UnitOrNonUnit diag, 
        AbstractDistMatrix<F>& APre,
  const AbstractDistMatrix<F>& UPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( APre.Height() != APre.Width() )
          LogicError("A must be square");
      if( UPre.Height() != UPre.Width() )
          LogicError("Triangular matrices must be square");
      if( APre.Height() != UPre.Height() )
          LogicError("A and U must be the same size");
    )
    const Int n = APre.Height();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    auto& A = AProx.Get();
    auto& U = UProx.GetLocked();

    // Temporary distributions
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g), U11_STAR_STAR(g),
                            X11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g), U12_STAR_VR(g);
    DistMatrix<F,STAR,MC  > U12_STAR_MC(g);
    DistMatrix<F,MC,  STAR> Z12Adj_MC_STAR(g);
    DistMatrix<F,MR,  STAR> U12Adj_MR_STAR(g), Z12Adj_MR_STAR(g);
    DistMatrix<F,VC,  STAR> U12Adj_VC_STAR(g);
    DistMatrix<F,MR,  MC  > Z12Adj_MR_MC(g);
    DistMatrix<F> Z12Adj(g), Y12(g);

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );
        auto A22 = A( ind2, ind2 );

        auto U11 = U( ind1, ind1 );
        auto U12 = U( ind1, ind2 );
        auto U22 = U( ind2, ind2 );

        // Y12 := U12 A22
        U12Adj_MR_STAR.AlignWith( A22 );
        Adjoint( U12, U12Adj_MR_STAR );
        U12Adj_VC_STAR.AlignWith( A22 );
        U12Adj_VC_STAR = U12Adj_MR_STAR;
        U12_STAR_MC.AlignWith( A22 );
        U12Adj_VC_STAR.AdjointPartialColAllGather( U12_STAR_MC );
        Z12Adj_MC_STAR.AlignWith( A22 );
        Z12Adj_MR_STAR.AlignWith( A22 );
        Z12Adj_MC_STAR.Resize( A12.Width(), nb );
        Z12Adj_MR_STAR.Resize( A12.Width(), nb );
        Zero( Z12Adj_MC_STAR );
        Zero( Z12Adj_MR_STAR );
        symm::LocalAccumulateRU
        ( ADJOINT, 
          F(1), A22, U12_STAR_MC, U12Adj_MR_STAR, 
          Z12Adj_MC_STAR, Z12Adj_MR_STAR );
        Z12Adj.AlignWith( A12 );
        Contract( Z12Adj_MC_STAR, Z12Adj );
        Z12Adj_MR_MC.AlignWith( A12 );
        Z12Adj_MR_MC = Z12Adj;
        AxpyContract( F(1), Z12Adj_MR_STAR, Z12Adj_MR_MC );
        Y12.AlignWith( A12 );
        Y12.Resize( nb, A12.Width() );
        Adjoint( Z12Adj_MR_MC.LockedMatrix(), Y12.Matrix() );

        // A12 := U11 A12
        A12_STAR_VR.AlignWith( A12 );
        A12_STAR_VR = A12;
        U11_STAR_STAR = U11;
        LocalTrmm
        ( LEFT, UPPER, NORMAL, diag, F(1), U11_STAR_STAR, A12_STAR_VR );
        A12 = A12_STAR_VR;

        // A12 := A12 + 1/2 Y12
        Axpy( F(1)/F(2), Y12, A12 );

        // A11 := U11 A11 U11'
        A11_STAR_STAR = A11;
        TwoSidedTrmm( UPPER, diag, A11_STAR_STAR, U11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A11 := A11 + (U12 A12' + A12 U12')
        A12_STAR_VR = A12;
        U12_STAR_VR.AlignWith( A12 );
        U12_STAR_VR = U12;
        X11_STAR_STAR.Resize( nb, nb );
        Zero( X11_STAR_STAR );
        Her2k
        ( UPPER, NORMAL,
          F(1), A12_STAR_VR.Matrix(), U12_STAR_VR.Matrix(),
          F(0), X11_STAR_STAR.Matrix() );
        AxpyContract( F(1), X11_STAR_STAR, A11 );

        // A12 := A12 + 1/2 Y12
        Axpy( F(1)/F(2), Y12, A12 );

        // A12 := A12 U22'
        Trmm( RIGHT, UPPER, ADJOINT, diag, F(1), U22, A12 );
    }
}

} // namespace twotrmm
} // namespace El

#endif // ifndef EL_TWOSIDEDTRMM_UVAR1_HPP
