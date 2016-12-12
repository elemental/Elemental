/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_INVERSE_HPD_CHOLESKYLVAR2_HPP
#define EL_INVERSE_HPD_CHOLESKYLVAR2_HPP

namespace El {
namespace hpd_inv {

// This approach is based upon a (conjugate)-transposition of the reordered
// Variant 2 algorithm from Fig. 9 in Bientinesi et al.'s "Families of
// Algorithms Related to the Inversion of a Symmetric Positive Definite Matrix".

template<typename Field>
void
CholeskyLVar2( Matrix<Field>& A )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("Nonsquare matrices cannot be triangular");
    )
    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A00 = A( ind0, ind0 );
        auto A10 = A( ind1, ind0 );
        auto A11 = A( ind1, ind1 );
        auto A20 = A( ind2, ind0 );
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );

        Cholesky( LOWER, A11 );
        Trsm( LEFT, LOWER, NORMAL, NON_UNIT, Field(1), A11, A10 );
        Trsm( RIGHT, LOWER, ADJOINT, NON_UNIT, Field(1), A11, A21 );
        Herk( LOWER, ADJOINT, Base<Field>(1), A10, Base<Field>(1), A00 );
        Gemm( NORMAL, NORMAL, Field(-1), A21, A10, Field(1), A20 );
        Herk( LOWER, NORMAL, Base<Field>(-1), A21, Base<Field>(1), A22 );
        Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, Field(1), A11, A10 );
        Trsm( RIGHT, LOWER, NORMAL, NON_UNIT, Field(-1), A11, A21 );
        TriangularInverse( LOWER, NON_UNIT, A11 );
        Trtrmm( LOWER, A11, true );
    }
}

template<typename Field>
void
CholeskyLVar2( AbstractDistMatrix<Field>& APre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( APre.Height() != APre.Width() )
          LogicError("Nonsquare matrices cannot be triangular");
    )

    DistMatrixReadWriteProxy<Field,Field,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Grid& g = A.Grid();
    DistMatrix<Field,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<Field,STAR,VR  > A10_STAR_VR(g);
    DistMatrix<Field,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<Field,STAR,MC  > A10_STAR_MC(g);
    DistMatrix<Field,STAR,MR  > A10_STAR_MR(g);
    DistMatrix<Field,STAR,MC  > A21Trans_STAR_MC(g);
    DistMatrix<Field,VR,  STAR> A21_VR_STAR(g);
    DistMatrix<Field,STAR,MR  > A21Adj_STAR_MR(g);

    const Int n = A.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A00 = A( ind0, ind0 );
        auto A10 = A( ind1, ind0 );
        auto A11 = A( ind1, ind1 );
        auto A20 = A( ind2, ind0 );
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );

        A21Trans_STAR_MC.AlignWith( A20 );
        A21_VR_STAR.AlignWith( A22 );
        A21Adj_STAR_MR.AlignWith( A22 );

        A11_STAR_STAR = A11;
        Cholesky( LOWER, A11_STAR_STAR );

        A10_STAR_VR.AlignWith( A00 );
        A10_STAR_VR = A10;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, NON_UNIT, Field(1), A11_STAR_STAR, A10_STAR_VR );

        A21_VC_STAR.AlignWith( A20 );
        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, NON_UNIT,
          Field(1), A11_STAR_STAR, A21_VC_STAR );

        A10_STAR_MC.AlignWith( A00 );
        A10_STAR_MC = A10_STAR_VR;
        A10_STAR_MR.AlignWith( A00 );
        A10_STAR_MR = A10_STAR_VR;
        LocalTrrk
        ( LOWER, ADJOINT,
          Field(1), A10_STAR_MC, A10_STAR_MR, Field(1), A00 );

        Transpose( A21_VC_STAR, A21Trans_STAR_MC );
        LocalGemm
        ( TRANSPOSE, NORMAL,
          Field(-1), A21Trans_STAR_MC, A10_STAR_MR, Field(1), A20 );

        A21_VR_STAR = A21_VC_STAR;
        Adjoint( A21_VR_STAR, A21Adj_STAR_MR );
        LocalTrrk
        ( LOWER, TRANSPOSE,
          Field(-1), A21Trans_STAR_MC, A21Adj_STAR_MR, Field(1), A22 );

        LocalTrsm
        ( LEFT, LOWER, ADJOINT, NON_UNIT,
          Field(1), A11_STAR_STAR, A10_STAR_VR );

        LocalTrsm
        ( RIGHT, LOWER, NORMAL, NON_UNIT,
          Field(-1), A11_STAR_STAR, A21_VC_STAR );

        LocalTriangularInverse( LOWER, NON_UNIT, A11_STAR_STAR );

        Trtrmm( LOWER, A11_STAR_STAR, true );

        A11 = A11_STAR_STAR;
        A10 = A10_STAR_VR;
        A21 = A21_VC_STAR;
    }
}

} // namespace hpd_inv
} // namespace El

#endif // ifndef EL_INVERSE_HPD_CHOLESKYLVAR2_HPP
