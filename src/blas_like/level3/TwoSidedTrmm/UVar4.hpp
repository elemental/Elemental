/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TWOSIDEDTRMM_UVAR4_HPP
#define EL_TWOSIDEDTRMM_UVAR4_HPP

namespace El {
namespace twotrmm {

// The only reason that a field is required is for the existence of 1/2, which
// is an artifact of this algorithm...
template<typename F> 
inline void
UVar4( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& U )
{
    DEBUG_ONLY(
      CSE cse("twotrmm::UVar4");
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

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A00 = A( ind0, ind0 );
        auto A01 = A( ind0, ind1 );
        auto A02 = A( ind0, ind2 );
        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        // Y01 := U01 A11
        Zeros( Y01, k, nb );
        Hemm( RIGHT, UPPER, F(1), A11, U01, F(0), Y01 );

        // A01 := A01 + 1/2 Y01
        Axpy( F(1)/F(2), Y01, A01 );

        // A00 := A00 + (U01 A01' + A01 U01')
        Her2k( UPPER, NORMAL, F(1), U01, A01, Base<F>(1), A00 );

        // A01 := A01 + 1/2 Y01
        Axpy( F(1)/F(2), Y01, A01 );

        // A01 := A01 U11'
        Trmm( RIGHT, UPPER, ADJOINT, diag, F(1), U11, A01 );

        // A02 := A02 + U01 A12
        Gemm( NORMAL, NORMAL, F(1), U01, A12, F(1), A02 );

        // A11 := U11 A11 U11'
        twotrmm::UUnb( diag, A11, U11 );

        // A12 := U11 A12
        Trmm( LEFT, UPPER, NORMAL, diag, F(1), U11, A12 );
    }
}

template<typename F> 
inline void
UVar4
( UnitOrNonUnit diag, 
        ElementalMatrix<F>& APre,
  const ElementalMatrix<F>& UPre )
{
    DEBUG_ONLY(
      CSE cse("twotrmm::UVar4");
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
    DistMatrix<F,STAR,MC  > A01Adj_STAR_MC(g), U01Adj_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A01Adj_STAR_MR(g), U01Adj_STAR_MR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,MR,  STAR> A12Adj_MR_STAR(g);
    DistMatrix<F,VC,  STAR> A01_VC_STAR(g), U01_VC_STAR(g), Y01_VC_STAR(g);
    DistMatrix<F,VR,  STAR> A01_VR_STAR(g), U01_VR_STAR(g);

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A00 = A( ind0, ind0 );
        auto A01 = A( ind0, ind1 );
        auto A02 = A( ind0, ind2 );
        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        // Y01 := U01 A11
        A11_STAR_STAR = A11;
        U01_VC_STAR.AlignWith( A00 );
        U01_VC_STAR = U01;
        Y01_VC_STAR.AlignWith( A01 );
        Zeros( Y01_VC_STAR, k, nb );
        Hemm
        ( RIGHT, UPPER, 
          F(1), A11_STAR_STAR.LockedMatrix(), U01_VC_STAR.LockedMatrix(), 
          F(0), Y01_VC_STAR.Matrix() );

        // A01 := A01 + 1/2 Y01
        A01_VC_STAR.AlignWith( A00 );
        A01_VC_STAR = A01;
        Axpy( F(1)/F(2), Y01_VC_STAR, A01_VC_STAR );

        // A00 := A00 + (U01 A01' + A01 U01')
        A01Adj_STAR_MC.AlignWith( A00 );
        U01Adj_STAR_MC.AlignWith( A00 );
        Adjoint( A01_VC_STAR, A01Adj_STAR_MC );
        Adjoint( U01_VC_STAR, U01Adj_STAR_MC );
        A01_VR_STAR.AlignWith( A00 );
        A01_VR_STAR = A01_VC_STAR;
        U01_VR_STAR.AlignWith( A00 );
        U01_VR_STAR = U01_VC_STAR;
        A01Adj_STAR_MR.AlignWith( A00 );
        U01Adj_STAR_MR.AlignWith( A00 );
        Adjoint( A01_VR_STAR, A01Adj_STAR_MR );
        Adjoint( U01_VR_STAR, U01Adj_STAR_MR );
        LocalTrr2k
        ( UPPER, ADJOINT, NORMAL, ADJOINT, NORMAL,
          F(1), U01Adj_STAR_MC, A01Adj_STAR_MR, 
          F(1), A01Adj_STAR_MC, U01Adj_STAR_MR, F(1), A00 );

        // A01 := A01 + 1/2 Y01
        Axpy( F(1)/F(2), Y01_VC_STAR, A01_VC_STAR );

        // A01 := A01 U11'
        U11_STAR_STAR = U11;
        LocalTrmm
        ( RIGHT, UPPER, ADJOINT, diag, F(1), U11_STAR_STAR, A01_VC_STAR );
        A01 = A01_VC_STAR;

        // A02 := A02 + U01 A12
        A12Adj_MR_STAR.AlignWith( A02 );
        Adjoint( A12, A12Adj_MR_STAR );
        LocalGemm
        ( ADJOINT, ADJOINT, F(1), U01Adj_STAR_MC, A12Adj_MR_STAR, F(1), A02 );

        // A11 := U11 A11 U11'
        TwoSidedTrmm( UPPER, diag, A11_STAR_STAR, U11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A12 := U11 A12
        Adjoint( A12Adj_MR_STAR, A12_STAR_VR );
        LocalTrmm
        ( LEFT, UPPER, NORMAL, diag, F(1), U11_STAR_STAR, A12_STAR_VR );
        A12 = A12_STAR_VR;
    }
}

} // namespace twotrmm
} // namespace El

#endif // ifndef EL_TWOSIDEDTRMM_UVAR4_HPP
