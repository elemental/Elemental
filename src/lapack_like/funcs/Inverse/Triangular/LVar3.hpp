/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_INVERSE_TRIANGULAR_LVAR3_HPP
#define EL_INVERSE_TRIANGULAR_LVAR3_HPP

namespace El {
namespace triang_inv {

template<typename Field>
void
LVar3Unb( UnitOrNonUnit diag, Matrix<Field>& L )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( L.Height() != L.Width() )
          LogicError("Nonsquare matrices cannot be triangular");
    )
    const Int n = L.Height();
    const Int ldl = L.LDim();
    Field* LBuffer = L.Buffer();
    for( Int j=0; j<n; ++j )
    {
        const Field lambda = ( diag==NON_UNIT ? LBuffer[j+j*ldl] : Field(1) );
        for( Int k=0; k<j; ++k )
            LBuffer[j+k*ldl] /= -lambda;
        blas::Geru
        ( n-(j+1), j, Field(1),
          &LBuffer[(j+1)+j*ldl], 1, &LBuffer[j], ldl,
          &LBuffer[j+1], ldl );
        if( diag == NON_UNIT )
        {
            for( Int k=j+1; k<n; ++k )
                LBuffer[k+j*ldl] /= lambda;
            LBuffer[j+j*ldl] = Field(1) / LBuffer[j+j*ldl];
        }
    }
}

template<typename Field>
void
LVar3( UnitOrNonUnit diag, Matrix<Field>& L )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( L.Height() != L.Width() )
          LogicError("Nonsquare matrices cannot be triangular");
    )
    const Int n = L.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );
        auto L20 = L( ind2, ind0 );
        auto L21 = L( ind2, ind1 );

        Trsm( LEFT, LOWER, NORMAL, diag, Field(-1), L11, L10 );
        Gemm( NORMAL, NORMAL, Field(1), L21, L10, Field(1), L20 );
        Trsm( RIGHT, LOWER, NORMAL, diag, Field(1), L11, L21 );
        LVar3Unb( diag, L11 );
    }
}

template<typename Field>
void
LVar3( UnitOrNonUnit diag, AbstractDistMatrix<Field>& LPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( LPre.Height() != LPre.Width() )
          LogicError("Nonsquare matrices cannot be triangular");
    )

    DistMatrixReadWriteProxy<Field,Field,MC,MR> LProx( LPre );
    auto& L = LProx.Get();

    const Grid& g = L.Grid();
    DistMatrix<Field,STAR,MR  > L10_STAR_MR(g);
    DistMatrix<Field,STAR,VR  > L10_STAR_VR(g);
    DistMatrix<Field,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<Field,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<Field,VC,  STAR> L21_VC_STAR(g);

    const Int n = L.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );
        auto L20 = L( ind2, ind0 );
        auto L21 = L( ind2, ind1 );

        L10_STAR_VR = L10;
        L11_STAR_STAR = L11;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, diag, Field(-1), L11_STAR_STAR, L10_STAR_VR );

        L21_MC_STAR.AlignWith( L20 );
        L21_MC_STAR = L21;
        L10_STAR_MR.AlignWith( L20 );
        L10_STAR_MR = L10_STAR_VR;
        LocalGemm
        ( NORMAL, NORMAL, Field(1), L21_MC_STAR, L10_STAR_MR, Field(1), L20 );
        L10 = L10_STAR_MR;

        L21_VC_STAR = L21_MC_STAR;
        LocalTrsm
        ( RIGHT, LOWER, NORMAL, diag, Field(1), L11_STAR_STAR, L21_VC_STAR );
        LocalTriangularInverse( LOWER, diag, L11_STAR_STAR );
        L11 = L11_STAR_STAR;
        L21 = L21_VC_STAR;
    }
}

} // namespace triang_inv
} // namespace El

#endif // ifndef EL_INVERSE_TRIANGULAR_LVAR3_HPP
