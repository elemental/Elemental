/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HESSENBERG_LOWER_BLOCKED_HPP
#define EL_HESSENBERG_LOWER_BLOCKED_HPP

#include "./LowerPanel.hpp"

namespace El {
namespace hessenberg {

template<typename F>
void LowerBlocked( Matrix<F>& A, Matrix<F>& householderScalars )
{
    EL_DEBUG_CSE
    const Int n = A.Height();
    householderScalars.Resize( Max(n-1,0), 1 );

    Matrix<F> UB1, V01, VB1, G11;

    const Int bsize = Blocksize();
    for( Int k=0; k<n-1; k+=bsize )
    {
        const Int nb = Min(bsize,n-1-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         indB( k,    n    ), indR( k, n ),
                         ind2( k+nb, n    );

        auto ABR = A( indB, indR );
        auto A22 = A( ind2, ind2 );

        auto householderScalars1 = householderScalars( ind1, ALL );
        UB1.Resize( n-k, nb );
        VB1.Resize( n-k, nb );
        G11.Resize( nb,  nb );
        hessenberg::LowerPanel( ABR, householderScalars1, UB1, VB1, G11 );

        auto AB0 = A( indB, ind0 );
        auto A2R = A( ind2, indR );
        auto U21 = UB1( IR(nb,END), ALL );
        auto V21 = VB1( IR(nb,END), ALL );

        // AB0 := AB0 - (UB1 inv(G11)^H UB1^H AB0)
        //      = AB0 - (UB1 ((AB0^H UB1) inv(G11))^H)
        // -------------------------------------------
        Gemm( ADJOINT, NORMAL, F(1), AB0, UB1, V01 );
        Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), G11, V01 );
        Gemm( NORMAL, ADJOINT, F(-1), UB1, V01, F(1), AB0 );

        // A2R := (A2R - U21 inv(G11)^H VB1^H)(I - UB1 inv(G11) UB1^H)
        // -----------------------------------------------------------
        // A2R := A2R - U21 inv(G11)^H VB1^H
        // (note: VB1 is overwritten)
        Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), G11, VB1 );
        Gemm( NORMAL, ADJOINT, F(-1), U21, VB1, F(1), A2R );
        // A2R := A2R - ((A2R UB1) inv(G11)) UB1^H
        Gemm( NORMAL, NORMAL, F(1), A2R, UB1, F(0), V21 );
        Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), G11, V21 );
        Gemm( NORMAL, ADJOINT, F(-1), V21, UB1, F(1), A2R );
    }
}

template<typename F>
void LowerBlocked
( AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<F>& householderScalarsPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(AssertSameGrids( APre, householderScalarsPre ))

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,STAR,STAR>
      householderScalarsProx( householderScalarsPre );
    auto& A = AProx.Get();
    auto& householderScalars = householderScalarsProx.Get();

    const Grid& g = A.Grid();
    const Int n = A.Height();
    householderScalars.Resize( Max(n-1,0), 1 );

    DistMatrix<F,MC,STAR> UB1_MC_STAR(g), V21_MC_STAR(g);
    DistMatrix<F,MR,STAR> V01_MR_STAR(g), VB1_MR_STAR(g), UB1_MR_STAR(g);
    DistMatrix<F,STAR,STAR> G11_STAR_STAR(g);

    const Int bsize = Blocksize();
    for( Int k=0; k<n-1; k+=bsize )
    {
        const Int nb = Min(bsize,n-1-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         indB( k,    n    ), indR( k, n ),
                         ind2( k+nb, n    );

        auto ABR = A( indB, indR );
        auto A22 = A( ind2, ind2 );

        auto householderScalars1 = householderScalars( ind1, ALL );
        UB1_MC_STAR.AlignWith( ABR );
        UB1_MR_STAR.AlignWith( ABR );
        VB1_MR_STAR.AlignWith( ABR );
        UB1_MC_STAR.Resize( n-k, nb );
        UB1_MR_STAR.Resize( n-k, nb );
        VB1_MR_STAR.Resize( n-k, nb );
        G11_STAR_STAR.Resize( nb, nb );
        hessenberg::LowerPanel
        ( ABR, householderScalars1, UB1_MC_STAR, UB1_MR_STAR, VB1_MR_STAR,
          G11_STAR_STAR );

        auto AB0 = A( indB, ind0 );
        auto A2R = A( ind2, indR );

        auto U21_MC_STAR = UB1_MC_STAR( IR(nb,END), ALL );

        // AB0 := AB0 - (UB1 inv(G11)^H UB1^H AB0)
        //      = AB0 - (UB1 ((AB0^H UB1) inv(G11))^H)
        // -------------------------------------------
        V01_MR_STAR.AlignWith( AB0 );
        Zeros( V01_MR_STAR, k, nb );
        LocalGemm( ADJOINT, NORMAL, F(1), AB0, UB1_MC_STAR, F(0), V01_MR_STAR );
        El::AllReduce( V01_MR_STAR, AB0.ColComm() );
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), G11_STAR_STAR, V01_MR_STAR );
        LocalGemm
        ( NORMAL, ADJOINT, F(-1), UB1_MC_STAR, V01_MR_STAR, F(1), AB0 );

        // A2R := (A2R - U21 inv(G11)^H VB1^H)(I - UB1 inv(G11) UB1^H)
        // -----------------------------------------------------------
        // A2R := A2R - U21 inv(G11)^H VB1^H
        // (note: VB1 is overwritten)
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), G11_STAR_STAR, VB1_MR_STAR );
        LocalGemm
        ( NORMAL, ADJOINT, F(-1), U21_MC_STAR, VB1_MR_STAR, F(1), A2R );
        // A2R := A2R - ((A2R UB1) inv(G11)) UB1^H
        V21_MC_STAR.AlignWith( A2R );
        Zeros( V21_MC_STAR, A2R.Height(), nb );
        LocalGemm( NORMAL, NORMAL, F(1), A2R, UB1_MR_STAR, F(0), V21_MC_STAR );
        El::AllReduce( V21_MC_STAR, A2R.RowComm() );
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), G11_STAR_STAR, V21_MC_STAR );
        LocalGemm
        ( NORMAL, ADJOINT, F(-1), V21_MC_STAR, UB1_MR_STAR, F(1), A2R );
    }
}

} // namespace hessenberg
} // namespace El

#endif // ifndef EL_HESSENBERG_LOWER_BLOCKED_HPP
