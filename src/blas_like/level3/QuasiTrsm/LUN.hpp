/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace quasitrsm {

// Left Upper Normal (Non)Unit Trsm
//   X := triu(U)^-1  X, or
//   X := triuu(U)^-1 X

template<typename F>
void LUNUnb( const Matrix<F>& U, Matrix<F>& X, bool checkIfSingular )
{
    EL_DEBUG_CSE
    const Int m = X.Height();
    const Int n = X.Width();
    typedef Base<F> Real;
    
    const F* UBuf = U.LockedBuffer();
          F* XBuf = X.Buffer();
    const Int ldu = U.LDim();
    const Int ldx = X.LDim();
    
    Int k=m-1;
    while( k >= 0 )
    {
        const bool in2x2 = ( k>0 && UBuf[k+(k-1)*ldu] != F(0) );
        if( in2x2 )
        {
            --k;
            // Solve the 2x2 linear systems via a 2x2 QR decomposition produced
            // by the Givens rotation
            //    | c        s | | U(k,  k) | = | gamma11 | 
            //    | -conj(s) c | | U(k+1,k) |   | 0       |
            //
            // and by also forming the right two entries of the 2x2 resulting
            // upper-triangular matrix, say gamma12 and gamma22
            //
            // Extract the 2x2 diagonal block, D
            const F delta11 = UBuf[ k   + k   *ldu];
            const F delta12 = UBuf[ k   +(k+1)*ldu];
            const F delta21 = UBuf[(k+1)+ k   *ldu];
            const F delta22 = UBuf[(k+1)+(k+1)*ldu];
            // Decompose D = Q R
            Real c; F s;
            const F gamma11 = Givens( delta11, delta21, c, s );
            const F gamma12 =        c*delta12 + s*delta22;
            const F gamma22 = -Conj(s)*delta12 + c*delta22;
            if( checkIfSingular )
            {
                // TODO: Instead check if values are too small in magnitude
                if( gamma11 == F(0) || gamma22 == F(0) )
                    LogicError("Singular diagonal block detected");
            }
            for( Int j=0; j<n; ++j )
            {
                F* xBuf = &XBuf[j*ldx];
                // Solve against Q
                const F chi1 = xBuf[k  ];
                const F chi2 = xBuf[k+1];
                xBuf[k  ] =        c*chi1 + s*chi2;
                xBuf[k+1] = -Conj(s)*chi1 + c*chi2;

                // Solve against R
                xBuf[k+1] /= gamma22;
                xBuf[k  ] -= gamma12*xBuf[k+1];
                xBuf[k  ] /= gamma11;

                // Update x0 := x0 - U01 x1
                blas::Axpy( k, -xBuf[k  ], &UBuf[ k   *ldu], 1, xBuf, 1 );
                blas::Axpy( k, -xBuf[k+1], &UBuf[(k+1)*ldu], 1, xBuf, 1 );
            }
        }
        else
        {
            if( checkIfSingular )
                if( UBuf[k+k*ldu] == F(0) )
                    LogicError("Singular diagonal entry detected");
            for( Int j=0; j<n; ++j )
            {
                F* xBuf = &XBuf[j*ldx];
                // Solve the 1x1 linear system
                xBuf[k] /= UBuf[k+k*ldu];

                // Update x0 := x0 - u01 chi_1
                blas::Axpy( k, -xBuf[k], &UBuf[k*ldu], 1, xBuf, 1 );
            }
        }
        --k;
    }
}

template<typename F>
void LUN( const Matrix<F>& U, Matrix<F>& X, bool checkIfSingular )
{
    EL_DEBUG_CSE
    const Int m = X.Height();
    const Int bsize = Blocksize();

    const Int kLast = LastOffset( m, bsize );
    Int k=kLast, kOld=m;
    while( true )
    {
        const bool in2x2 = ( k>0 && U.Get(k,k-1) != F(0) );
        if( in2x2 )
            --k;
        const Int nb = kOld-k;

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        auto X0 = X( ind0, ALL );
        auto X1 = X( ind1, ALL );

        LUNUnb( U11, X1, checkIfSingular );
        Gemm( NORMAL, NORMAL, F(-1), U01, X1, F(1), X0 );

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
}

template<typename F>
void LUNLarge
( const AbstractDistMatrix<F>& UPre,
        AbstractDistMatrix<F>& XPre, 
  bool checkIfSingular )
{
    EL_DEBUG_CSE
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<F,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    const Int kLast = LastOffset( m, bsize );
    Int k=kLast, kOld=m;
    while( true )
    {
        const bool in2x2 = ( k>0 && U.Get(k,k-1) != F(0) );
        if( in2x2 )
            --k;
        const Int nb = kOld-k;

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        auto X0 = X( ind0, ALL );
        auto X1 = X( ind1, ALL );

        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        X1_STAR_VR    = X1;  // X1[* ,VR] <- X1[MC,MR]
        
        // X1[* ,VR] := U11^-1[* ,* ] X1[* ,VR]
        LocalQuasiTrsm
        ( LEFT, UPPER, NORMAL, F(1), U11_STAR_STAR, X1_STAR_VR,
          checkIfSingular );

        X1_STAR_MR.AlignWith( X0 );
        X1_STAR_MR  = X1_STAR_VR; // X1[* ,MR]  <- X1[* ,VR]
        X1          = X1_STAR_MR; // X1[MC,MR] <- X1[* ,MR]
        U01_MC_STAR.AlignWith( X0 );
        U01_MC_STAR = U01;        // U01[MC,* ] <- U01[MC,MR]

        // X0[MC,MR] -= U01[MC,* ] X1[* ,MR]
        LocalGemm( NORMAL, NORMAL, F(-1), U01_MC_STAR, X1_STAR_MR, F(1), X0 );

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
}

template<typename F>
void LUNMedium
( const AbstractDistMatrix<F>& UPre, AbstractDistMatrix<F>& XPre, 
  bool checkIfSingular )
{
    EL_DEBUG_CSE
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<F,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,MR,  STAR> X1Trans_MR_STAR(g);

    const Int kLast = LastOffset( m, bsize );
    Int k=kLast, kOld=m;
    while( true )
    {
        const bool in2x2 = ( k>0 && U.Get(k,k-1) != F(0) );
        if( in2x2 )
            --k;
        const Int nb = kOld-k;

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        auto X0 = X( ind0, ALL );
        auto X1 = X( ind1, ALL );

        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        X1Trans_MR_STAR.AlignWith( X0 );
        Transpose( X1, X1Trans_MR_STAR );
        
        // X1^T[MR,* ] := X1^T[MR,* ] U11^-T[* ,* ]
        //              = (U11^-1[* ,* ] X1[* ,MR])^T
        LocalQuasiTrsm
        ( RIGHT, UPPER, TRANSPOSE,
          F(1), U11_STAR_STAR, X1Trans_MR_STAR, checkIfSingular );
        Transpose( X1Trans_MR_STAR, X1 );

        U01_MC_STAR.AlignWith( X0 );
        U01_MC_STAR = U01;  // U01[MC,* ] <- U01[MC,MR]

        // X0[MC,MR] -= U01[MC,* ] X1[* ,MR]
        LocalGemm
        ( NORMAL, TRANSPOSE, F(-1), U01_MC_STAR, X1Trans_MR_STAR, F(1), X0 );

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
}

template<typename F,Dist colDist>
void LUNSmall
( const DistMatrix<F,colDist,STAR>& U,
        DistMatrix<F,colDist,STAR>& X,
  bool checkIfSingular )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( U, X );
      if( U.Height() != U.Width() || U.Width() != X.Height() )
          LogicError
          ("Nonconformal: \n",DimsString(U,"U"),"\n",DimsString(X,"X"));
      if( U.ColAlign() != X.ColAlign() )
          LogicError("U and X are assumed to be aligned");
    )
    const Int m = X.Height();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g), X1_STAR_STAR(g);

    const Int kLast = LastOffset( m, bsize );
    Int k=kLast, kOld=m;
    while( true )
    {
        const bool in2x2 = ( k>0 && U.Get(k,k-1) != F(0) );
        if( in2x2 )
            --k;
        const Int nb = kOld-k;

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        auto X0 = X( ind0, ALL );
        auto X1 = X( ind1, ALL );

        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[VC,* ]
        X1_STAR_STAR = X1;   // X1[* ,* ] <- X1[VC,* ]
        
        // X1[* ,* ] := U11^-1[* ,* ] X1[* ,* ]
        LocalQuasiTrsm
        ( LEFT, UPPER, NORMAL, 
          F(1), U11_STAR_STAR, X1_STAR_STAR, checkIfSingular );
        X1 = X1_STAR_STAR;

        // X0[VC,* ] -= U01[VC,* ] X1[* ,* ]
        LocalGemm( NORMAL, NORMAL, F(-1), U01, X1_STAR_STAR, F(1), X0 );

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
}

} // namespace quasitrsm
} // namespace El
