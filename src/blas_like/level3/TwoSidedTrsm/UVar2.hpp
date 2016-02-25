/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TWOSIDEDTRSM_UVAR2_HPP
#define EL_TWOSIDEDTRSM_UVAR2_HPP

namespace El {
namespace twotrsm {

template<typename F> 
inline void
UVar2( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& U )
{
    DEBUG_ONLY(
      CSE cse("twotrsm::UVar2");
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

        // Y01 := A00 U01
        Zeros( Y01, A01.Height(), A01.Width() );
        Hemm( LEFT, UPPER, F(1), A00, U01, F(0), Y01 );

        // A01 := A01 - 1/2 Y01
        Axpy( F(-1)/F(2), Y01, A01 );
        
        // A11 := A11 - (U01' A01 + A01' U01)
        Her2k( UPPER, ADJOINT, F(-1), U01, A01, F(1), A11 );

        // A11 := inv(U11)' A11 inv(U11)
        twotrsm::UUnb( diag, A11, U11 );

        // A12 := A12 - A02' U01
        Gemm( ADJOINT, NORMAL, F(-1), A02, U01, F(1), A12 );

        // A12 := inv(U11)' A12
        Trsm( LEFT, UPPER, ADJOINT, diag, F(1), U11, A12 );
        
        // A01 := A01 - 1/2 Y01
        Axpy( F(-1)/F(2), Y01, A01 );

        // A01 := A01 inv(U11)
        Trsm( RIGHT, UPPER, NORMAL, diag, F(1), U11, A01 );
    }
}

// This routine has only partially been optimized. The ReduceScatter operations
// need to be (conjugate-)transposed in order to play nice with cache.
template<typename F> 
inline void
UVar2
( UnitOrNonUnit diag, 
        ElementalMatrix<F>& APre,
  const ElementalMatrix<F>& UPre )
{
    DEBUG_ONLY(
      CSE cse("twotrsm::UVar2");
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
    DistMatrix<F,MC,  STAR> A01_MC_STAR(g), F01_MC_STAR(g), U01_MC_STAR(g);
    DistMatrix<F,VC,  STAR> A01_VC_STAR(g);
    DistMatrix<F,VR,  STAR> U01_VR_STAR(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g), U11_STAR_STAR(g);
    DistMatrix<F,STAR,MR  > U01Adj_STAR_MR(g), X11_STAR_MR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,MR,  STAR> X12Adj_MR_STAR(g), Y01_MR_STAR(g);
    DistMatrix<F,MR,  MC  > X12Adj_MR_MC(g), Y01_MR_MC(g);
    DistMatrix<F> X11(g), Y01(g);

    Matrix<F> X12Local;

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

        // Y01 := A00 U01
        U01_MC_STAR.AlignWith( A00 );
        U01_MC_STAR = U01;
        U01_VR_STAR.AlignWith( A00 );
        U01_VR_STAR = U01_MC_STAR;
        U01Adj_STAR_MR.AlignWith( A00 );
        U01_VR_STAR.AdjointPartialColAllGather( U01Adj_STAR_MR );
        Y01_MR_STAR.AlignWith( A00 );
        F01_MC_STAR.AlignWith( A00 );
        Zeros( Y01_MR_STAR, k, nb );
        Zeros( F01_MC_STAR, k, nb );
        symm::LocalAccumulateLU
        ( ADJOINT, 
          F(1), A00, U01_MC_STAR, U01Adj_STAR_MR, F01_MC_STAR, Y01_MR_STAR );
        Contract( Y01_MR_STAR, Y01_MR_MC );
        Y01.AlignWith( A01 );
        Y01 = Y01_MR_MC;
        AxpyContract( F(1), F01_MC_STAR, Y01 );

        // X11 := U01' A01
        X11_STAR_MR.AlignWith( U01 );
        LocalGemm( ADJOINT, NORMAL, F(1), U01_MC_STAR, A01, X11_STAR_MR );

        // A01 := A01 - Y01
        A01 -= Y01;
        A01_MC_STAR.AlignWith( U01 );
        A01_MC_STAR = A01;
        
        // A11 := A11 - triu(X11 + A01' U01) = A11 - (U01 A01 + A01' U01)
        LocalGemm( ADJOINT, NORMAL, F(1), A01_MC_STAR, U01, F(1), X11_STAR_MR );
        X11.AlignWith( A11 );
        Contract( X11_STAR_MR, X11 );
        MakeTrapezoidal( UPPER, X11 );
        A11 -= X11;

        // A01 := A01 inv(U11)
        U11_STAR_STAR = U11;
        A01_VC_STAR = A01_MC_STAR;
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, diag, F(1), U11_STAR_STAR, A01_VC_STAR );
        A01 = A01_VC_STAR;

        // A11 := inv(U11)' A11 inv(U11)
        A11_STAR_STAR = A11;
        TwoSidedTrsm( UPPER, diag, A11_STAR_STAR, U11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A12 := A12 - A02' U01
        X12Adj_MR_STAR.AlignWith( A02 );
        LocalGemm( ADJOINT, NORMAL, F(1), A02, U01_MC_STAR, X12Adj_MR_STAR );
        X12Adj_MR_MC.AlignWith( A12 );
        Contract( X12Adj_MR_STAR, X12Adj_MR_MC );
        Adjoint( X12Adj_MR_MC.LockedMatrix(), X12Local );
        A12.Matrix() -= X12Local;

        // A12 := inv(U11)' A12
        A12_STAR_VR = A12;
        LocalTrsm
        ( LEFT, UPPER, ADJOINT, diag, F(1), U11_STAR_STAR, A12_STAR_VR );
        A12 = A12_STAR_VR;
    }
}

} // namespace twotrsm
} // namespace El

#endif // ifndef EL_TWOSIDEDTRSM_UVAR2_HPP
