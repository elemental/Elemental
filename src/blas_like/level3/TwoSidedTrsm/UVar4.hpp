/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TWOSIDEDTRSM_UVAR4_HPP
#define EL_TWOSIDEDTRSM_UVAR4_HPP

namespace El {
namespace twotrsm {

template<typename F> 
void UVar4( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& U )
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
    Matrix<F> Y12;

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
        auto A22 = A( ind2, ind2 );

        auto U11 = U( ind1, ind1 );
        auto U12 = U( ind1, ind2 );

        // A01 := A01 inv(U11)
        Trsm( RIGHT, UPPER, NORMAL, diag, F(1), U11, A01 );

        // A11 := inv(U11)' A11 inv(U11)
        twotrsm::UUnb( diag, A11, U11 );

        // A02 := A02 - A01 U12
        Gemm( NORMAL, NORMAL, F(-1), A01, U12, F(1), A02 );

        // Y12 := A11 U12
        Y12.Resize( A12.Height(), A12.Width() );
        Zero( Y12 );
        Hemm( LEFT, UPPER, F(1), A11, U12, F(0), Y12 );

        // A12 := inv(U11)' A12
        Trsm( LEFT, UPPER, ADJOINT, diag, F(1), U11, A12 );

        // A12 := A12 - 1/2 Y12
        Axpy( F(-1)/F(2), Y12, A12 );

        // A22 := A22 - (A12' U12 + U12' A12)
        Her2k( UPPER, ADJOINT, F(-1), A12, U12, Base<F>(1), A22 );

        // A12 := A12 - 1/2 Y12
        Axpy( F(-1)/F(2), Y12, A12 );
    }
}

template<typename F> 
void UVar4
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
    DistMatrix<F,STAR,MC  > A01Trans_STAR_MC(g), A12_STAR_MC(g), U12_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A12_STAR_MR(g);
    DistMatrix<F,STAR,VC  > A12_STAR_VC(g), U12_STAR_VC(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g), U12_STAR_VR(g), Y12_STAR_VR(g);
    DistMatrix<F,MR,  STAR> U12Trans_MR_STAR(g);
    DistMatrix<F,VC,  STAR> A01_VC_STAR(g);
    DistMatrix<F,VR,  STAR> U12Trans_VR_STAR(g);

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
        auto A22 = A( ind2, ind2 );

        auto U11 = U( ind1, ind1 );
        auto U12 = U( ind1, ind2 );

        // A01 := A01 inv(U11)
        A01_VC_STAR.AlignWith( A02 );
        A01_VC_STAR = A01;
        U11_STAR_STAR = U11;
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, diag, F(1), U11_STAR_STAR, A01_VC_STAR );
        A01 = A01_VC_STAR;

        // A11 := inv(U11)' A11 inv(U11)
        A11_STAR_STAR = A11;
        TwoSidedTrsm( UPPER, diag, A11_STAR_STAR, U11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A02 := A02 - A01 U12
        A01Trans_STAR_MC.AlignWith( A02 );
        Transpose( A01_VC_STAR, A01Trans_STAR_MC );
        U12Trans_MR_STAR.AlignWith( A02 );
        Transpose( U12, U12Trans_MR_STAR );
        LocalGemm
        ( TRANSPOSE, TRANSPOSE, 
          F(-1), A01Trans_STAR_MC, U12Trans_MR_STAR, F(1), A02 );

        // Y12 := A11 U12
        U12Trans_VR_STAR.AlignWith( A02 );
        U12Trans_VR_STAR = U12Trans_MR_STAR;
        U12_STAR_VR.AlignWith( A02 );
        U12_STAR_VR.Resize( nb, A12.Width() );
        Zero( U12_STAR_VR ); 
        Transpose( U12Trans_VR_STAR.Matrix(), U12_STAR_VR.Matrix() );
        Y12_STAR_VR.AlignWith( A12 );
        Y12_STAR_VR.Resize( nb, A12.Width() );
        Zero( Y12_STAR_VR );
        Hemm
        ( LEFT, UPPER, 
          F(1), A11_STAR_STAR.Matrix(), U12_STAR_VR.Matrix(), 
          F(0), Y12_STAR_VR.Matrix() );

        // A12 := inv(U11)' A12
        A12_STAR_VR.AlignWith( A22 );
        A12_STAR_VR = A12;
        LocalTrsm
        ( LEFT, UPPER, ADJOINT, diag, F(1), U11_STAR_STAR, A12_STAR_VR );

        // A12 := A12 - 1/2 Y12
        Axpy( F(-1)/F(2), Y12_STAR_VR, A12_STAR_VR );

        // A22 := A22 - (A12' U12 + U12' A12)
        A12_STAR_MR.AlignWith( A22 );
        A12_STAR_MR = A12_STAR_VR;
        A12_STAR_VC.AlignWith( A22 );
        A12_STAR_VC = A12_STAR_VR;
        U12_STAR_VC.AlignWith( A22 );
        U12_STAR_VC = U12_STAR_VR;
        A12_STAR_MC.AlignWith( A22 );
        A12_STAR_MC = A12_STAR_VC;
        U12_STAR_MC.AlignWith( A22 );
        U12_STAR_MC = U12_STAR_VC;
        LocalTrr2k
        ( UPPER, ADJOINT, TRANSPOSE, ADJOINT, NORMAL,
          F(-1), A12_STAR_MC, U12Trans_MR_STAR,
          F(-1), U12_STAR_MC, A12_STAR_MR,
          F(1), A22 );

        // A12 := A12 - 1/2 Y12
        Axpy( F(-1)/F(2), Y12_STAR_VR, A12_STAR_VR );
        A12 = A12_STAR_VR;
    }
}

} // namespace twotrsm
} // namespace El

#endif // ifndef EL_TWOSIDEDTRSM_UVAR4_HPP
