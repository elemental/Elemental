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

template<typename F>
void
LVar3Unb( UnitOrNonUnit diag, Matrix<F>& L )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( L.Height() != L.Width() )
          LogicError("Nonsquare matrices cannot be triangular");
    )
    const Int n = L.Height();
    const Int ldl = L.LDim();
    F* LBuffer = L.Buffer();
    for( Int j=0; j<n; ++j )
    {
        const F lambda = ( diag==NON_UNIT ? LBuffer[j+j*ldl] : F(1) );
        for( Int k=0; k<j; ++k )
            LBuffer[j+k*ldl] /= -lambda;
        blas::Geru
        ( n-(j+1), j, F(1),
          &LBuffer[(j+1)+j*ldl], 1, &LBuffer[j], ldl, 
          &LBuffer[j+1], ldl );
        if( diag == NON_UNIT )
        {
            for( Int k=j+1; k<n; ++k )
                LBuffer[k+j*ldl] /= lambda;
            LBuffer[j+j*ldl] = F(1) / LBuffer[j+j*ldl];
        }
    }
}

template<typename F>
void
LVar3( UnitOrNonUnit diag, Matrix<F>& L )
{
    DEBUG_CSE
    DEBUG_ONLY(
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

        Trsm( LEFT, LOWER, NORMAL, diag, F(-1), L11, L10 );
        Gemm( NORMAL, NORMAL, F(1), L21, L10, F(1), L20 );
        Trsm( RIGHT, LOWER, NORMAL, diag, F(1), L11, L21 );
        LVar3Unb( diag, L11 );
    }
}

template<typename F>
void
LVar3( UnitOrNonUnit diag, ElementalMatrix<F>& LPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( LPre.Height() != LPre.Width() )
          LogicError("Nonsquare matrices cannot be triangular");
    )

    DistMatrixReadWriteProxy<F,F,MC,MR> LProx( LPre );
    auto& L = LProx.Get();

    const Grid& g = L.Grid();
    DistMatrix<F,STAR,MR  > L10_STAR_MR(g);
    DistMatrix<F,STAR,VR  > L10_STAR_VR(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<F,VC,  STAR> L21_VC_STAR(g);

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
        ( LEFT, LOWER, NORMAL, diag, F(-1), L11_STAR_STAR, L10_STAR_VR );

        L21_MC_STAR.AlignWith( L20 );
        L21_MC_STAR = L21;
        L10_STAR_MR.AlignWith( L20 );
        L10_STAR_MR = L10_STAR_VR;
        LocalGemm
        ( NORMAL, NORMAL, F(1), L21_MC_STAR, L10_STAR_MR, F(1), L20 );
        L10 = L10_STAR_MR;

        L21_VC_STAR = L21_MC_STAR;
        LocalTrsm
        ( RIGHT, LOWER, NORMAL, diag, F(1), L11_STAR_STAR, L21_VC_STAR );
        LocalTriangularInverse( LOWER, diag, L11_STAR_STAR );
        L11 = L11_STAR_STAR;
        L21 = L21_VC_STAR;
    }
}

} // namespace triang_inv
} // namespace El

#endif // ifndef EL_INVERSE_TRIANGULAR_LVAR3_HPP
