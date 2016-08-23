/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRSTRM_LLN_HPP
#define EL_TRSTRM_LLN_HPP

namespace El {
namespace trstrm {

template<typename F>
void LLNUnb( UnitOrNonUnit diag, F alpha, const Matrix<F>& L, Matrix<F>& X )
{
    DEBUG_CSE
    const bool isUnit = ( diag==UNIT );
    const Int n = L.Height();
    const Int LLDim = L.LDim();
    const Int XLDim = X.LDim();
    const F* LBuffer = L.LockedBuffer();
    F* XBuffer = X.Buffer();

    // X := alpha X
    if( alpha != F(1) )
        for( Int j=0; j<n; ++j ) 
            for( Int i=j; i<n; ++i )
                XBuffer[i+j*XLDim] *= alpha;

    for( Int i=0; i<n; ++i )
    {
        if( !isUnit )
        {
            const F lambda11 = LBuffer[i+i*LLDim];
            for( Int j=0; j<i; ++j )
                XBuffer[i+j*XLDim] /= lambda11;
            XBuffer[i+i*XLDim] /= lambda11;
        }

        const Int l21Height = n - (i+1);
        const F* l21 = &LBuffer[(i+1)+i*LLDim];
        const F* x1L = &XBuffer[i];
        F* X2L = &XBuffer[i+1];
        blas::Geru( l21Height, i+1, F(-1), l21, 1, x1L, XLDim, X2L, XLDim );
    }
}

template<typename F>
void LLN
( UnitOrNonUnit diag,
  F alpha,
  const Matrix<F>& L,
        Matrix<F>& X,
  bool checkIfSingular=true )
{
    DEBUG_CSE
    const Int n = L.Height();
    const Int bsize = Blocksize();

    Matrix<F> Z11;

    ScaleTrapezoid( alpha, LOWER, X );
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto L11 = L( ind1, ind1 );
        auto L21 = L( ind2, ind1 );

        auto X10 = X( ind1, ind0 );
        auto X11 = X( ind1, ind1 );
        auto X20 = X( ind2, ind0 );
        auto X21 = X( ind2, ind1 );

        Trsm( LEFT, LOWER, NORMAL, diag, F(1), L11, X10, checkIfSingular );
        trstrm::LLNUnb( diag, F(1), L11, X11 );
        Gemm( NORMAL, NORMAL, F(-1), L21, X10, F(1), X20 );
        Z11 = X11;
        MakeTrapezoidal( LOWER, Z11 );
        Gemm( NORMAL, NORMAL, F(-1), L21, Z11, F(1), X21 );
    }
}

template<typename F>
void LLN
( UnitOrNonUnit diag, 
  F alpha,
  const AbstractDistMatrix<F>& LPre,
        AbstractDistMatrix<F>& XPre,
  bool checkIfSingular )
{
    DEBUG_CSE
    const Int n = LPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> LProx( LPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& L = LProx.GetLocked();
    auto& X = XProx.Get();

    // Temporary distributions
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g), X11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<F,STAR,MR  > X10_STAR_MR(g), X11_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X10_STAR_VR(g);

    ScaleTrapezoid( alpha, LOWER, X );
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto L11 = L( ind1, ind1 );
        auto L21 = L( ind2, ind1 );

        auto X10 = X( ind1, ind0 );
        auto X11 = X( ind1, ind1 );
        auto X20 = X( ind2, ind0 );
        auto X21 = X( ind2, ind1 );

        L11_STAR_STAR = L11; 
        X11_STAR_STAR = X11; 
        X10_STAR_VR = X10;

        LocalTrsm
        ( LEFT, LOWER, NORMAL, diag, F(1), L11_STAR_STAR, X10_STAR_VR,
          checkIfSingular );
        Trstrm
        ( LEFT, LOWER, NORMAL, diag, F(1), L11_STAR_STAR, X11_STAR_STAR,
          checkIfSingular );
        X11 = X11_STAR_STAR;
        X11_STAR_MR.AlignWith( X21 );
        X11_STAR_MR = X11_STAR_STAR;
        MakeTrapezoidal( LOWER, X11_STAR_MR );

        X10_STAR_MR.AlignWith( X20 );
        X10_STAR_MR = X10_STAR_VR;
        X10 = X10_STAR_MR;
        L21_MC_STAR.AlignWith( X20 );
        L21_MC_STAR = L21;
        
        LocalGemm
        ( NORMAL, NORMAL, F(-1), L21_MC_STAR, X10_STAR_MR, F(1), X20 );
        LocalGemm
        ( NORMAL, NORMAL, F(-1), L21_MC_STAR, X11_STAR_MR, F(1), X21 );
    }
}

} // namespace trstrm
} // namespace El

#endif // ifndef EL_TRSTRM_LLN_HPP
