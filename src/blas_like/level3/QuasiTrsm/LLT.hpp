/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace quasitrsm {

// Left Lower (Conjugate)Transpose (Non)Unit QuasiTrsm
//   X := tril(L)^-T,
//   X := tril(L)^-H,
//   X := trilu(L)^-T, or
//   X := trilu(L)^-H

template<typename F>
void LLTUnb
( bool conjugate, const Matrix<F>& L, Matrix<F>& X, bool checkIfSingular )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Int m = X.Height();
    const Int n = X.Width();

    const F* LBuf = L.LockedBuffer();
          F* XBuf = X.Buffer();
    const Int ldl = L.LDim();
    const Int ldx = X.LDim();

    if( conjugate )
        Conjugate( X );

    Int k=m-1;
    while( k >= 0 )
    {
        const bool in2x2 = ( k>0 && LBuf[(k-1)+k*ldl] != F(0) );
        if( in2x2 )
        {
            --k;
            // Solve the 2x2 linear systems via a 2x2 LQ decomposition produced
            // by the Givens rotation
            //    | L(k,k) L(k,k+1) | | c -conj(s) | = | gamma11 0 |
            //                        | s    c     |
            // and by also forming the bottom two entries of the 2x2 resulting
            // lower-triangular matrix, say gamma21 and gamma22
            //
            // Extract the 2x2 diagonal block, D
            const F delta11 = LBuf[   k +   k *ldl];
            const F delta12 = LBuf[   k +(k+1)*ldl];
            const F delta21 = LBuf[(k+1)+   k *ldl];
            const F delta22 = LBuf[(k+1)+(k+1)*ldl];
            // Decompose D = L Q
            Real c; F s;
            const F gamma11 = Givens( delta11, delta12, c, s );
            const F gamma21 =        c*delta21 + s*delta22;
            const F gamma22 = -Conj(s)*delta21 + c*delta22;
            if( checkIfSingular )
            {
                // TODO: Instead check if values are too small in magnitude
                if( gamma11 == F(0) || gamma22 == F(0) )
                    LogicError("Singular diagonal block detected");
            }
            for( Int j=0; j<n; ++j )
            {
                F* xBuf = &XBuf[j*ldx];
                // Solve against Q^T
                const F chi1 = xBuf[k  ];
                const F chi2 = xBuf[k+1];
                xBuf[k  ] =        c*chi1 + s*chi2;
                xBuf[k+1] = -Conj(s)*chi1 + c*chi2;

                // Solve against R^T
                xBuf[k+1] /= gamma22;
                xBuf[k  ] -= gamma21*xBuf[k+1];
                xBuf[k  ] /= gamma11;

                // Update x0 := x0 - L10^T x1
                blas::Axpy( k, -xBuf[k  ], &LBuf[k  ], ldl, xBuf, 1 );
                blas::Axpy( k, -xBuf[k+1], &LBuf[k+1], ldl, xBuf, 1 );
            }
        }
        else
        {
            if( checkIfSingular )
                if( LBuf[k+k*ldl] == F(0) )
                    LogicError("Singular diagonal entry detected");
            for( Int j=0; j<n; ++j )
            {
                F* xBuf = &XBuf[j*ldx];
                // Solve the 1x1 linear system
                xBuf[k] /= LBuf[k+k*ldl];

                // Update x0 := x0 - l10^T chi_1
                blas::Axpy( k, -xBuf[k], &LBuf[k], ldl, xBuf, 1 );
            }
        }
        --k;
    }
    if( conjugate )
        Conjugate( X );
}

template<typename F>
void LLT
( Orientation orientation, 
  const Matrix<F>& L, Matrix<F>& X, bool checkIfSingular )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( orientation == NORMAL )
          LogicError("Expected (Conjugate)Transpose option");
    )
    const Int m = X.Height();
    const Int bsize = Blocksize();

    const bool conjugate = ( orientation==ADJOINT );

    const Int kLast = LastOffset( m, bsize );
    Int k=kLast, kOld=m;
    while( true )
    {
        const bool in2x2 = ( k>0 && L.Get(k-1,k) != F(0) );
        if( in2x2 )
            --k;
        const Int nb = kOld-k;

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );

        auto X0 = X( ind0, ALL );
        auto X1 = X( ind1, ALL );

        LLTUnb( conjugate, L11, X1, checkIfSingular );
        Gemm( orientation, NORMAL, F(-1), L10, X1, F(1), X0 );

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
}

// width(X) >> p
template<typename F>
void LLTLarge
( Orientation orientation, 
  const AbstractDistMatrix<F>& LPre,
        AbstractDistMatrix<F>& XPre,
  bool checkIfSingular )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( orientation == NORMAL )
          LogicError("Expected (Conjugate)Transpose option");
    )
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> LProx( LPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& L = LProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<F,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    const Int kLast = LastOffset( m, bsize );
    Int k=kLast, kOld=m;
    while( true )
    {
        const bool in2x2 = ( k>0 && L.Get(k-1,k) != F(0) );
        if( in2x2 )
            --k;
        const Int nb = kOld-k;

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );

        auto X0 = X( ind0, ALL );
        auto X1 = X( ind1, ALL );

        // X1[* ,VR] := L11^-[T/H][* ,* ] X1[* ,VR]
        L11_STAR_STAR = L11; 
        X1_STAR_VR    = X1;
        LocalQuasiTrsm
        ( LEFT, LOWER, orientation, F(1), L11_STAR_STAR, X1_STAR_VR,
          checkIfSingular );

        X1_STAR_MR.AlignWith( X0 );
        X1_STAR_MR  = X1_STAR_VR; // X1[* ,MR] <- X1[* ,VR]
        X1          = X1_STAR_MR; // X1[MC,MR] <- X1[* ,MR]
        L10_STAR_MC.AlignWith( X0 );
        L10_STAR_MC = L10;        // L10[* ,MC] <- L10[MC,MR]

        // X0[MC,MR] -= (L10[* ,MC])^(T/H) X1[* ,MR]
        //            = L10^[T/H][MC,* ] X1[* ,MR]
        LocalGemm
        ( orientation, NORMAL, F(-1), L10_STAR_MC, X1_STAR_MR, F(1), X0 );

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
}

// width(X) ~= p
template<typename F>
void LLTMedium
( Orientation orientation, 
  const AbstractDistMatrix<F>& LPre,
        AbstractDistMatrix<F>& XPre,
  bool checkIfSingular )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( orientation == NORMAL )
          LogicError("Expected (Conjugate)Transpose option");
    )
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> LProx( LPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& L = LProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<F,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,MR,  STAR> X1Trans_MR_STAR(g);

    const Int kLast = LastOffset( m, bsize );
    Int k=kLast, kOld=m;
    while( true )
    {
        const bool in2x2 = ( k>0 && L.Get(k-1,k) != F(0) );
        if( in2x2 )
            --k;
        const Int nb = kOld-k;

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );

        auto X0 = X( ind0, ALL );
        auto X1 = X( ind1, ALL );

        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[MC,MR]
        // X1[* ,MR] <- X1[MC,MR]
        X1Trans_MR_STAR.AlignWith( X0 );
        Transpose( X1, X1Trans_MR_STAR, (orientation==ADJOINT) );

        // X1[* ,MR] := L11^-[T/H][* ,* ] X1[* ,MR]
        // X1^[T/H][MR,* ] := X1^[T/H][MR,* ] L11^-1[* ,* ]
        LocalQuasiTrsm
        ( RIGHT, LOWER, NORMAL,
          F(1), L11_STAR_STAR, X1Trans_MR_STAR, checkIfSingular );

        Transpose( X1Trans_MR_STAR, X1, (orientation==ADJOINT) );
        L10_STAR_MC.AlignWith( X0 );
        L10_STAR_MC = L10; // L10[* ,MC] <- L10[MC,MR]

        // X0[MC,MR] -= (L10[* ,MC])^[T/H] X1[* ,MR]
        //            = L10^[T/H][MC,* ] X1[* ,MR]
        LocalGemm
        ( orientation, orientation, 
          F(-1), L10_STAR_MC, X1Trans_MR_STAR, F(1), X0 );

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
}

// width(X) << p
template<typename F,Dist colDist>
void LLTSmall
( Orientation orientation, 
  const DistMatrix<F,colDist,STAR>& L,
        DistMatrix<F,colDist,STAR>& X,
  bool checkIfSingular )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( L, X );
      if( orientation == NORMAL )
          LogicError("Expected (Conjugate)Transpose option");
      if( L.Height() != L.Width() || L.Height() != X.Height() )
          LogicError
          ("Nonconformal: \n",DimsString(L,"L"),"\n",DimsString(X,"X"));
      if( L.ColAlign() != X.ColAlign() )
          LogicError("L and X must be aligned");
    )
    const Int m = X.Height();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> Z1_STAR_STAR(g);

    const Int kLast = LastOffset( m, bsize );
    Int k=kLast, kOld=m;
    while( true )
    {
        const bool in2x2 = ( k>0 && L.Get(k-1,k) != F(0) );
        if( in2x2 )
            --k;
        const Int nb = kOld-k;

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, m    );

        auto L11 = L( ind1, ind1 );
        auto L21 = L( ind2, ind1 );

        auto X1 = X( ind1, ALL );
        auto X2 = X( ind2, ALL );

        // X1 -= L21' X2
        LocalGemm( orientation, NORMAL, F(-1), L21, X2, Z1_STAR_STAR );
        axpy::util::UpdateWithLocalData( F(1), X1, Z1_STAR_STAR );
        El::AllReduce( Z1_STAR_STAR, X1.DistComm() );

        // X1 := L11^-1 X1
        L11_STAR_STAR = L11;
        LocalQuasiTrsm
        ( LEFT, LOWER, orientation, F(1), L11_STAR_STAR, Z1_STAR_STAR,
          checkIfSingular );
        X1 = Z1_STAR_STAR;

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
}

template<typename F,Dist rowDist>
void LLTSmall
( Orientation orientation, 
  const DistMatrix<F,STAR,rowDist>& L,
        DistMatrix<F,rowDist,STAR>& X,
  bool checkIfSingular )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( L, X );
      if( orientation == NORMAL )
          LogicError("Expected (Conjugate)Transpose option");
      if( L.Height() != L.Width() || L.Height() != X.Height() )
          LogicError
          ("Nonconformal: \n",DimsString(L,"L"),"\n",DimsString(X,"X"));
      if( L.RowAlign() != X.ColAlign() )
          LogicError("L and X must be aligned");
    )
    const Int m = X.Height();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g), X1_STAR_STAR(g);

    const Int kLast = LastOffset( m, bsize );
    Int k=kLast, kOld=m;
    while( true )
    {
        const bool in2x2 = ( k>0 && L.Get(k-1,k) != F(0) );
        if( in2x2 )
            --k;
        const Int nb = kOld-k;

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );

        auto X0 = X( ind0, ALL );
        auto X1 = X( ind1, ALL );

        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[* ,VR]
        X1_STAR_STAR = X1;   // X1[* ,* ] <- X1[VR,* ]

        // X1[* ,* ] := L11^-[T/H][* ,* ] X1[* ,* ]
        LocalQuasiTrsm
        ( LEFT, LOWER, orientation,
          F(1), L11_STAR_STAR, X1_STAR_STAR, checkIfSingular );

        X1 = X1_STAR_STAR;

        // X0[VR,* ] -= L10[* ,VR]^(T/H) X1[* ,* ]
        LocalGemm( orientation, NORMAL, F(-1), L10, X1_STAR_STAR, F(1), X0 );

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
}

} // namespace quasitrsm
} // namespace El
