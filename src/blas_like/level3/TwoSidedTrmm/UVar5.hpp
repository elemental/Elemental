/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TWOSIDEDTRMM_UVAR5_HPP
#define EL_TWOSIDEDTRMM_UVAR5_HPP

namespace El {
namespace twotrmm {

// The only requirement that this is a field comes from the necessity for 
// the existence of 1/2, which is artifact of the algorithm...
template<typename F> 
void UVar5( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& U )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
    Matrix<F> Y01;

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto A00 = A( ind0, ind0 );
        auto A01 = A( ind0, ind1 );
        auto A11 = A( ind1, ind1 );

        auto U00 = U( ind0, ind0 );
        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        // Y01 := U01 A11
        Y01.Resize( k, nb );
        Zero( Y01 );
        Hemm( RIGHT, UPPER, F(1), A11, U01, F(0), Y01 );

        // A01 := U00 A01
        Trmm( LEFT, UPPER, NORMAL, diag, F(1), U00, A01 );

        // A01 := A01 + 1/2 Y01
        Axpy( F(1)/F(2), Y01, A01 );

        // A00 := A00 + (U01 A01' + A01 U01')
        Her2k( UPPER, NORMAL, F(1), U01, A01, Base<F>(1), A00 );

        // A01 := A01 + 1/2 Y01
        Axpy( F(1)/F(2), Y01, A01 );

        // A01 := A01 U11'
        Trmm( RIGHT, UPPER, ADJOINT, diag, F(1), U11, A01 );

        // A11 := U11 A11 U11'
        twotrmm::UUnb( diag, A11, U11 );
    }
}

template<typename F> 
void UVar5
( UnitOrNonUnit diag, 
        AbstractDistMatrix<F>& APre,
  const AbstractDistMatrix<F>& UPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g), U11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> A01_MC_STAR(g), U01_MC_STAR(g);
    DistMatrix<F,MR,  STAR> A01_MR_STAR(g), U01_MR_STAR(g);
    DistMatrix<F,VC,  STAR> A01_VC_STAR(g), U01_VC_STAR(g), Y01_VC_STAR(g);
    DistMatrix<F> Y01(g);

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto A00 = A( ind0, ind0 );
        auto A01 = A( ind0, ind1 );
        auto A11 = A( ind1, ind1 );

        auto U00 = U( ind0, ind0 );
        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        // Y01 := U01 A11
        A11_STAR_STAR = A11;
        U01_VC_STAR.AlignWith( A00 );
        U01_VC_STAR = U01;
        Y01_VC_STAR.AlignWith( A01 );
        Y01_VC_STAR.Resize( k, nb );
        Zero( Y01_VC_STAR );
        Hemm
        ( RIGHT, UPPER,
          F(1), A11_STAR_STAR.Matrix(), U01_VC_STAR.Matrix(),
          F(0), Y01_VC_STAR.Matrix() );
        Y01.AlignWith( A01 );
        Y01 = Y01_VC_STAR;

        // A01 := U00 A01
        Trmm( LEFT, UPPER, NORMAL, diag, F(1), U00, A01 );

        // A01 := A01 + 1/2 Y01
        Axpy( F(1)/F(2), Y01, A01 );

        // A00 := A00 + (U01 A01' + A01 U01')
        A01_MC_STAR.AlignWith( A00 );
        A01_MC_STAR = A01;
        U01_MC_STAR.AlignWith( A00 );
        U01_MC_STAR = U01;
        A01_VC_STAR.AlignWith( A00 );
        A01_VC_STAR = A01_MC_STAR;
        A01_MR_STAR.AlignWith( A00 );
        A01_MR_STAR = A01_VC_STAR;
        U01_MR_STAR.AlignWith( A00 );
        U01_MR_STAR = U01_MC_STAR;
        LocalTrr2k
        ( UPPER, NORMAL, ADJOINT, NORMAL, ADJOINT,
          F(1), U01_MC_STAR, A01_MR_STAR, 
          F(1), A01_MC_STAR, U01_MR_STAR,
          F(1), A00 );

        // A01 := A01 + 1/2 Y01
        Axpy( F(1)/F(2), Y01_VC_STAR, A01_VC_STAR );

        // A01 := A01 U11'
        U11_STAR_STAR = U11;
        LocalTrmm
        ( RIGHT, UPPER, ADJOINT, diag, F(1), U11_STAR_STAR, A01_VC_STAR );
        A01 = A01_VC_STAR;

        // A11 := U11 A11 U11'
        TwoSidedTrmm( UPPER, diag, A11_STAR_STAR, U11_STAR_STAR );
        A11 = A11_STAR_STAR;
    }
}

} // namespace twotrmm
} // namespace El

#endif // ifndef EL_TWOSIDEDTRMM_UVAR5_HPP
