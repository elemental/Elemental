/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace msquasitrsm {

template<typename F>
void LLNUnb( const Matrix<F>& L, const Matrix<F>& shifts, Matrix<F>& X )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Int m = X.Height();
    const Int n = X.Width();

    const F* LBuf = L.LockedBuffer();
          F* XBuf = X.Buffer();
    const Int ldl = L.LDim();
    const Int ldx = X.LDim();

    Int k=0;
    while( k < m )
    {
        const bool in2x2 = ( k+1<m && LBuf[k+(k+1)*ldl] != F(0) );
        if( in2x2 ) 
        {
            // Solve the 2x2 linear systems via 2x2 LQ decompositions produced
            // by the Givens rotation
            //    | L(k,k)-shift L(k,k+1) | | c -conj(s) | = | gamma11 0 |
            //                              | s    c     |
            // and by also forming the bottom two entries of the 2x2 resulting
            // lower-triangular matrix, say gamma21 and gamma22
            //
            // Extract the constant part of the 2x2 diagonal block, D
            const F delta12 = LBuf[ k   +(k+1)*ldl];
            const F delta21 = LBuf[(k+1)+ k   *ldl];
            for( Int j=0; j<n; ++j )
            {
                const F delta11 = LBuf[ k   + k   *ldl] - shifts.Get(j,0);
                const F delta22 = LBuf[(k+1)+(k+1)*ldl] - shifts.Get(j,0);
                // Decompose D = L Q
                Real c; F s;
                const F gamma11 = Givens( delta11, delta12, c, s );
                const F gamma21 =        c*delta21 + s*delta22;
                const F gamma22 = -Conj(s)*delta21 + c*delta22;

                F* xBuf = &XBuf[j*ldx];

                // Solve against L
                xBuf[k  ] /= gamma11;
                xBuf[k+1] -= gamma21*xBuf[k];
                xBuf[k+1] /= gamma22;

                // Solve against Q
                const F chi1 = xBuf[k  ];
                const F chi2 = xBuf[k+1];
                xBuf[k  ] = c*chi1 - Conj(s)*chi2;
                xBuf[k+1] = s*chi1 +       c*chi2;

                // Update x2 := x2 - L21 x1
                blas::Axpy
                ( m-(k+2), -xBuf[k  ], 
                  &LBuf[(k+2)+ k   *ldl], 1, &xBuf[k+2], 1 );
                blas::Axpy
                ( m-(k+2), -xBuf[k+1], 
                  &LBuf[(k+2)+(k+1)*ldl], 1, &xBuf[k+2], 1 );
            }

            k += 2;
        }
        else
        {
            for( Int j=0; j<n; ++j )
            {
                F* xBuf = &XBuf[j*ldx];
                xBuf[k] /= LBuf[k+k*ldl] - shifts.Get(j,0);
                blas::Axpy
                ( m-(k+1), -xBuf[k], &LBuf[(k+1)+k*ldl], 1, &xBuf[k+1], 1 );
            }
            k += 1;
        }
    }
}

template<typename F>
void LLN( const Matrix<F>& L, const Matrix<F>& shifts, Matrix<F>& X )
{
    EL_DEBUG_CSE
    const Int m = X.Height();
    const Int bsize = Blocksize();

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = ( k+nbProp<m && L.Get(k+nbProp-1,k+nbProp) != F(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, m    );

        auto L11 = L( ind1, ind1 );
        auto L21 = L( ind2, ind1 );

        auto X1 = X( ind1, ALL );
        auto X2 = X( ind2, ALL );

        LLNUnb( L11, shifts, X1 );
        Gemm( NORMAL, NORMAL, F(-1), L21, X1, F(1), X2 );
    }
}

// For large numbers of RHS's, e.g., width(X) >> p
template<typename F>
void LLNLarge
( const AbstractDistMatrix<F>& LPre,
  const AbstractDistMatrix<F>& shiftsPre, 
        AbstractDistMatrix<F>& XPre )
{
    EL_DEBUG_CSE
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> LProx( LPre );
    DistMatrixReadProxy<F,F,VR,STAR> shiftsProx( shiftsPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& L = LProx.GetLocked();
    auto& X = XProx.Get();
    auto& shifts = shiftsProx.GetLocked();

    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = ( k+nbProp<m && L.Get(k+nbProp-1,k+nbProp) != F(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, m    );

        auto L11 = L( ind1, ind1 );
        auto L21 = L( ind2, ind1 );

        auto X1 = X( ind1, ALL );
        auto X2 = X( ind2, ALL );

        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[MC,MR]
        X1_STAR_VR.AlignWith( shifts );
        X1_STAR_VR    = X1;  // X1[* ,VR] <- X1[MC,MR]

        // X1[* ,VR] := L11^-1[* ,* ] X1[* ,VR]
        LocalMultiShiftQuasiTrsm
        ( LEFT, LOWER, NORMAL, F(1), L11_STAR_STAR, shifts, X1_STAR_VR );

        X1_STAR_MR.AlignWith( X2 );
        X1_STAR_MR  = X1_STAR_VR; // X1[* ,MR]  <- X1[* ,VR]
        X1          = X1_STAR_MR; // X1[MC,MR] <- X1[* ,MR]
        L21_MC_STAR.AlignWith( X2 );
        L21_MC_STAR = L21;        // L21[MC,* ] <- L21[MC,MR]
        
        // X2[MC,MR] -= L21[MC,* ] X1[* ,MR]
        LocalGemm( NORMAL, NORMAL, F(-1), L21_MC_STAR, X1_STAR_MR, F(1), X2 );
    }
}

// For medium numbers of RHS's, e.g., width(X) ~= p
template<typename F>
void LLNMedium
( const AbstractDistMatrix<F>& LPre,
  const AbstractDistMatrix<F>& shiftsPre, 
        AbstractDistMatrix<F>& XPre )
{
    EL_DEBUG_CSE
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> LProx( LPre );
    DistMatrixReadProxy<F,F,VR,STAR> shiftsProx( shiftsPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& L = LProx.GetLocked();
    auto& X = XProx.Get();
    auto& shifts = shiftsProx.GetLocked();

    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<F,MR,  STAR> X1Trans_MR_STAR(g);

    DistMatrix<F,MR,  STAR> shifts_MR_STAR( shifts ),
                            shifts_MR_STAR_Align(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = ( k+nbProp<m && L.Get(k+nbProp-1,k+nbProp) != F(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, m    );

        auto L11 = L( ind1, ind1 );
        auto L21 = L( ind2, ind1 );

        auto X1 = X( ind1, ALL );
        auto X2 = X( ind2, ALL );

        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[MC,MR]
        X1Trans_MR_STAR.AlignWith( X2 );
        Transpose( X1, X1Trans_MR_STAR );

        // X1^T[MR,* ] := X1^T[MR,* ] L11^-T[* ,* ]
        //              = (L11^-1[* ,* ] X1[* ,MR])^T
        shifts_MR_STAR_Align.AlignWith( X1Trans_MR_STAR );
        shifts_MR_STAR_Align = shifts_MR_STAR; 
        LocalMultiShiftQuasiTrsm
        ( RIGHT, LOWER, TRANSPOSE,
          F(1), L11_STAR_STAR, shifts_MR_STAR_Align, X1Trans_MR_STAR );

        Transpose( X1Trans_MR_STAR, X1 );
        L21_MC_STAR.AlignWith( X2 );
        L21_MC_STAR = L21;                   // L21[MC,* ] <- L21[MC,MR]
        
        // X2[MC,MR] -= L21[MC,* ] X1[* ,MR]
        LocalGemm
        ( NORMAL, TRANSPOSE, F(-1), L21_MC_STAR, X1Trans_MR_STAR, F(1), X2 );
    }
}

// For small numbers of RHS's, e.g., width(X) < p
template<typename F,Dist colDist,Dist shiftColDist,Dist shiftRowDist>
void LLNSmall
( const DistMatrix<F,     colDist,STAR        >& L, 
  const DistMatrix<F,shiftColDist,shiftRowDist>& shifts, 
        DistMatrix<F,     colDist,STAR        >& X )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( L.ColAlign() != X.ColAlign() )
          LogicError("L and X are assumed to be aligned");
    )
    const Int m = X.Height();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g), X1_STAR_STAR(g),
                            shifts_STAR_STAR(shifts);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = ( k+nbProp<m && L.Get(k+nbProp-1,k+nbProp) != F(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, m    );

        auto L11 = L( ind1, ind1 );
        auto L21 = L( ind2, ind1 );

        auto X1 = X( ind1, ALL );
        auto X2 = X( ind2, ALL );

        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[VC,* ]
        X1_STAR_STAR = X1;   // X1[* ,* ] <- X1[VC,* ]

        // X1[* ,* ] := (L11[* ,* ])^-1 X1[* ,* ]
        LocalMultiShiftQuasiTrsm
        ( LEFT, LOWER, NORMAL, 
          F(1), L11_STAR_STAR, shifts_STAR_STAR, X1_STAR_STAR );

        // X2[VC,* ] -= L21[VC,* ] X1[* ,* ]
        LocalGemm( NORMAL, NORMAL, F(-1), L21, X1_STAR_STAR, F(1), X2 );
    }
}

} // namespace msquasitrsm
} // namespace El
