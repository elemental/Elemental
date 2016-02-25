/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TWOSIDEDTRSM_UVAR1_HPP
#define EL_TWOSIDEDTRSM_UVAR1_HPP

namespace El {
namespace twotrsm {

template<typename F> 
inline void
UVar1( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& U )
{
    DEBUG_ONLY(
      CSE cse("twotrsm::UVar1");
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

        // Y01 := A00 U01
        Zeros( Y01, A01.Height(), A01.Width() );
        Hemm( LEFT, UPPER, F(1), A00, U01, F(0), Y01 );

        // A01 := inv(U00)' A01
        Trsm( LEFT, UPPER, ADJOINT, diag, F(1), U00, A01 );

        // A01 := A01 - 1/2 Y01
        Axpy( F(-1)/F(2), Y01, A01 );

        // A11 := A11 - (U01' A01 + A01' U01)
        Her2k( UPPER, ADJOINT, F(-1), U01, A01, F(1), A11 );

        // A11 := inv(U11)' A11 inv(U11)
        twotrsm::UUnb( diag, A11, U11 );

        // A01 := A01 - 1/2 Y01
        Axpy( F(-1)/F(2), Y01, A01 );

        // A01 := A01 inv(U11)
        Trsm( RIGHT, UPPER, NORMAL, diag, F(1), U11, A01 );
    }
}

template<typename F> 
inline void
UVar1
( UnitOrNonUnit diag, 
        ElementalMatrix<F>& APre,
  const ElementalMatrix<F>& UPre )
{
    DEBUG_ONLY(
      CSE cse("twotrsm::UVar1");
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
    DistMatrix<F,STAR,MR  > U01Adj_STAR_MR(g);
    DistMatrix<F,VC,  STAR> A01_VC_STAR(g), U01_VC_STAR(g);
    DistMatrix<F,VR,  STAR> U01_VR_STAR(g);
    DistMatrix<F,MC,  STAR> U01_MC_STAR(g), Z01_MC_STAR(g);
    DistMatrix<F,MR,  STAR> Z01_MR_STAR(g);
    DistMatrix<F,MR,  MC  > Z01_MR_MC(g);
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

        // Y01 := A00 U01
        U01_MC_STAR.AlignWith( A00 );
        U01_MC_STAR = U01;
        U01_VR_STAR.AlignWith( A00 );
        U01_VR_STAR = U01_MC_STAR;
        U01Adj_STAR_MR.AlignWith( A00 );
        U01_VR_STAR.AdjointPartialColAllGather( U01Adj_STAR_MR );
        Z01_MC_STAR.AlignWith( A00 );
        Z01_MR_STAR.AlignWith( A00 );
        Zeros( Z01_MC_STAR, A01.Height(), A01.Width() );
        Zeros( Z01_MR_STAR, A01.Height(), A01.Width() );
        symm::LocalAccumulateLU
        ( ADJOINT, 
          F(1), A00, U01_MC_STAR, U01Adj_STAR_MR, Z01_MC_STAR, Z01_MR_STAR );
        Z01_MR_MC.AlignWith( A01 );
        Contract( Z01_MR_STAR, Z01_MR_MC );
        Y01.AlignWith( A01 );
        Y01 = Z01_MR_MC;
        AxpyContract( F(1), Z01_MC_STAR, Y01 );

        // A01 := inv(U00)' A01
        //
        // This is the bottleneck because A01 only has blocksize columns
        Trsm( LEFT, UPPER, ADJOINT, diag, F(1), U00, A01 );

        // A01 := A01 - 1/2 Y01
        Axpy( F(-1)/F(2), Y01, A01 );

        // A11 := A11 - (U01' A01 + A01' U01)
        A01_VC_STAR.AlignWith( A01 );
        A01_VC_STAR = A01;
        U01_VC_STAR.AlignWith( A00 );
        U01_VC_STAR = U01_MC_STAR;
        Zeros( X11_STAR_STAR, A11.Height(), A11.Width() );
        Her2k
        ( UPPER, ADJOINT,
          F(-1), A01_VC_STAR.Matrix(), U01_VC_STAR.Matrix(),
          F(0), X11_STAR_STAR.Matrix() );
        AxpyContract( F(1), X11_STAR_STAR, A11 );

        // A11 := inv(U11)' A11 inv(U11)
        A11_STAR_STAR = A11;
        U11_STAR_STAR = U11;
        TwoSidedTrsm( UPPER, diag, A11_STAR_STAR, U11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A01 := A01 - 1/2 Y01
        Axpy( F(-1)/F(2), Y01, A01 );

        // A01 := A01 inv(U11)
        A01_VC_STAR = A01;
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, diag, F(1), U11_STAR_STAR, A01_VC_STAR );
        A01 = A01_VC_STAR;
    }
}

} // namespace twotrsm
} // namespace El

#endif // ifndef EL_TWOSIDEDTRSM_UVAR1_HPP
