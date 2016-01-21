/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TWOSIDEDTRSM_UVAR3_HPP
#define EL_TWOSIDEDTRSM_UVAR3_HPP

namespace El {
namespace twotrsm {

template<typename F> 
inline void
UVar3( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& U )
{
    DEBUG_ONLY(
      CSE cse("twotrsm::UVar4");
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( U.Height() != U.Width() )
          LogicError("Triangular matrices must be square");
      if( A.Height() != U.Height() )
          LogicError("A and U must be the same size");
    )
    const Int n = A.Height();
    const Int bsize = Blocksize();

    // We will use an entire extra matrix as temporary storage. If this is not
    // acceptable, use UVar4 instead.
    Matrix<F> Y;
    Zeros( Y, n, n );

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A01 = A( ind0, ind1 );
        auto A02 = A( ind0, ind2 );
        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );
        auto U12 = U( ind1, ind2 );

        auto Y01 = Y( ind0, ind1 );
        auto Y02 = Y( ind0, ind2 );
        auto Y12 = Y( ind1, ind2 );
 
        // A01 := A01 - 1/2 Y01
        Axpy( F(-1)/F(2), Y01, A01 );

        // A11 := A11 - (A01' U01 + U01' A01)
        Her2k( UPPER, ADJOINT, F(-1), A01, U01, F(1), A11 );

        // A11 := inv(U11)' A11 inv(U11)
        twotrsm::UUnb( diag, A11, U11 );

        // A12 := A12 - U01' A02
        Gemm( ADJOINT, NORMAL, F(-1), U01, A02, F(1), A12 );

        // A12 := inv(U11)' A12
        Trsm( LEFT, UPPER, ADJOINT, diag, F(1), U11, A12 );

        // A01 := A01 - 1/2 Y01
        Axpy( F(-1)/F(2), Y01, A01 );

        // A01 := A01 inv(U11)
        Trsm( RIGHT, UPPER, NORMAL, diag, F(1), U11, A01 );

        // Y02 := Y02 + A01 U12
        Gemm( NORMAL, NORMAL, F(1), A01, U12, F(1), Y02 );

        // Y12 := Y12 + A11 U12
        Hemm( LEFT, UPPER, F(1), A11, U12, F(1), Y12 );

        // Y12 := Y12 + A01' U02
        Gemm( ADJOINT, NORMAL, F(1), A01, U02, F(1), Y12 );
    }
}

template<typename F> 
inline void
UVar3
( UnitOrNonUnit diag, 
        ElementalMatrix<F>& APre,
  const ElementalMatrix<F>& UPre )
{
    DEBUG_ONLY(
      CSE cse("twotrsm::UVar4");
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
    DistMatrix<F,MC,  STAR> A01_MC_STAR(g), U01_MC_STAR(g), A11_MC_STAR(g);
    DistMatrix<F,MR,  STAR> U12Adj_MR_STAR(g);
    DistMatrix<F,VC,  STAR> A01_VC_STAR(g), U01_VC_STAR(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g), U11_STAR_STAR(g), 
                            X11_STAR_STAR(g);
    DistMatrix<F,STAR,MR  > X12_STAR_MR(g), Z12_STAR_MR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);

    // We will use an entire extra matrix as temporary storage. If this is not
    // acceptable, use UVar4 instead.
    DistMatrix<F> Y(g);
    Y.AlignWith( A );
    Zeros( Y, n, n );

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A01 = A( ind0, ind1 );
        auto A02 = A( ind0, ind2 );
        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );
        auto U12 = U( ind1, ind2 );

        auto Y01 = Y( ind0, ind1 );
        auto Y02 = Y( ind0, ind2 );
        auto Y12 = Y( ind1, ind2 );

        // A01 := A01 - 1/2 Y01
        Axpy( F(-1)/F(2), Y01, A01 );

        // A11 := A11 - (A01' U01 + U01' A01)
        A01_VC_STAR.AlignWith( A01 );
        U01_VC_STAR.AlignWith( A01 );
        A01_VC_STAR = A01;
        U01_VC_STAR = U01;
        Zeros( X11_STAR_STAR, A11.Height(), A11.Width() );
        Her2k
        ( UPPER, ADJOINT, 
          F(1), A01_VC_STAR.Matrix(), U01_VC_STAR.Matrix(),
          F(0), X11_STAR_STAR.Matrix() );
        MakeTrapezoidal( UPPER, X11_STAR_STAR );
        AxpyContract( F(-1), X11_STAR_STAR, A11 );

        // A11 := inv(U11)' A11 inv(U11)
        A11_STAR_STAR = A11;
        U11_STAR_STAR = U11;
        TwoSidedTrsm( UPPER, diag, A11_STAR_STAR, U11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A12 := A12 - U01' A02
        U01_MC_STAR.AlignWith( A01 );
        U01_MC_STAR = U01;
        X12_STAR_MR.AlignWith( A02 );
        LocalGemm( ADJOINT, NORMAL, F(1), U01_MC_STAR, A02, X12_STAR_MR );
        AxpyContract( F(-1), X12_STAR_MR, A12 );

        // A12 := inv(U11)' A12
        A12_STAR_VR.AlignWith( A12 );
        A12_STAR_VR = A12;
        LocalTrsm
        ( LEFT, UPPER, ADJOINT, diag, F(1), U11_STAR_STAR, A12_STAR_VR );
        A12 = A12_STAR_VR;

        // A01 := A01 - 1/2 Y01
        Axpy( F(-1)/F(2), Y01, A01 );

        // A01 := A01 inv(U11)
        A01_VC_STAR = A01;
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, diag, F(1), U11_STAR_STAR, A01_VC_STAR );
        A01 = A01_VC_STAR;

        // Y02 := Y02 + A01 U12
        A01_MC_STAR.AlignWith( A01 );
        A01_MC_STAR = A01;
        U12Adj_MR_STAR.AlignWith( Y12 );
        U12.AdjointColAllGather( U12Adj_MR_STAR );
        LocalGemm
        ( NORMAL, ADJOINT, F(1), A01_MC_STAR, U12Adj_MR_STAR, F(1), Y02 );

        // Y12 := Y12 + A11 U12
        MakeHermitian( UPPER, A11_STAR_STAR );
        A11_MC_STAR.AlignWith( Y12 );
        A11_MC_STAR = A11_STAR_STAR;
        LocalGemm
        ( NORMAL, ADJOINT, F(1), A11_MC_STAR, U12Adj_MR_STAR, F(0), Y12 );

        // Y12 := Y12 + A01' U02
        Z12_STAR_MR.AlignWith( U02 );
        LocalGemm( ADJOINT, NORMAL, F(1), A01_MC_STAR, U02, Z12_STAR_MR );
        AxpyContract( F(1), Z12_STAR_MR, Y12 );
    }
}

} // namespace twotrsm
} // namespace El

#endif // ifndef EL_TWOSIDEDTRSM_UVAR3_HPP
