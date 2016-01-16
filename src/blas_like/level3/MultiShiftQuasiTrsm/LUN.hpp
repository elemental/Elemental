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
inline void
LUNUnb( const Matrix<F>& U, const Matrix<F>& shifts, Matrix<F>& X )
{
    DEBUG_ONLY(
      CSE cse("msquasitrsm::LUNUnb");
      if( U.Height() != U.Width() )
          LogicError("U should be square");
      if( X.Height() != U.Height() )
          LogicError("X should be the same height as U's size");
    )

    const Int m = X.Height();
    const Int n = X.Width();
    typedef Base<F> Real;
    
    const F* UBuf = U.LockedBuffer();
          F* XBuf = X.Buffer();
    const Int ldu = U.LDim();
    const Int ldX = X.LDim();
    
    Int k=m-1;
    while( k >= 0 )
    {
        const bool in2x2 = ( k>0 && UBuf[k+(k-1)*ldu] != F(0) );
        if( in2x2 )
        {
            --k;
            // Solve the 2x2 linear systems via 2x2 QR decompositions produced
            // by the Givens rotation
            //    | c        s | | U(k,  k)-shift | = | gamma11 | 
            //    | -conj(s) c | | U(k+1,k)       |   | 0       |
            //
            // and by also forming the right two entries of the 2x2 resulting
            // upper-triangular matrix, say gamma12 and gamma22
            //
            // Extract the constant part of the 2x2 diagonal block, D
            const F delta12 = UBuf[ k   +(k+1)*ldu];
            const F delta21 = UBuf[(k+1)+ k   *ldu];
            for( Int j=0; j<n; ++j )
            {
                const F delta11 = UBuf[ k   + k   *ldu] - shifts.Get(j,0);
                const F delta22 = UBuf[(k+1)+(k+1)*ldu] - shifts.Get(j,0);
                // Decompose D = Q R
                Real c; F s;
                const F gamma11 = blas::Givens( delta11, delta21, &c, &s );
                const F gamma12 =        c*delta12 + s*delta22;
                const F gamma22 = -Conj(s)*delta12 + c*delta22;

                F* xBuf = &XBuf[j*ldX];
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
            for( Int j=0; j<n; ++j )
            {
                F* xBuf = &XBuf[j*ldX];
                // Solve the 1x1 linear system
                xBuf[k] /= UBuf[k+k*ldu] - shifts.Get(j,0);

                // Update x0 := x0 - u01 chi_1
                blas::Axpy( k, -xBuf[k], &UBuf[k*ldu], 1, xBuf, 1 );
            }
        }
        --k;
    }
}

template<typename Real>
inline void
LUNUnb
( const Matrix<Real>& U, 
  const Matrix<Complex<Real>>& shifts, 
        Matrix<Real>& XReal, Matrix<Real>& XImag )
{
    DEBUG_ONLY(
      CSE cse("msquasitrsm::LUNUnb");
      if( U.Height() != U.Width() )
          LogicError("U should be square");
      if( XReal.Height() != XImag.Height() ||
          XReal.Width()  != XImag.Width() )
          LogicError("XReal and XImag should be the same size");
      if( XReal.Height() != U.Height() )
          LogicError("X should be the same height as U's size");
    )
    const Int m = XReal.Height();
    const Int n = XReal.Width();
    typedef Complex<Real> C;
    
    const Real* UBuf = U.LockedBuffer();
          Real* XRealBuf = XReal.Buffer();
          Real* XImagBuf = XImag.Buffer();
    const Int ldu = U.LDim();
    const Int ldXReal = XReal.LDim();
    const Int ldXImag = XImag.LDim();
    
    Int k=m-1;
    while( k >= 0 )
    {
        const bool in2x2 = ( k>0 && UBuf[k+(k-1)*ldu] != Real(0) );
        if( in2x2 )
        {
            --k;
            // Solve the 2x2 linear systems via 2x2 QR decompositions produced
            // by the Givens rotation
            //    | c        s | | U(k,  k)-shift | = | gamma11 | 
            //    | -conj(s) c | | U(k+1,k)       |   | 0       |
            //
            // and by also forming the right two entries of the 2x2 resulting
            // upper-triangular matrix, say gamma12 and gamma22
            //
            // Extract the constant part of the 2x2 diagonal block, D
            const Real delta12 = UBuf[ k   +(k+1)*ldu];
            const Real delta21 = UBuf[(k+1)+ k   *ldu];
            for( Int j=0; j<n; ++j )
            {
                const C delta11 = UBuf[ k   + k   *ldu] - shifts.Get(j,0);
                const C delta22 = UBuf[(k+1)+(k+1)*ldu] - shifts.Get(j,0);
                // Decompose D = Q R
                Real c; C s;
                const C gamma11 = blas::Givens( delta11, C(delta21), &c, &s );
                const C gamma12 =        c*delta12 + s*delta22;
                const C gamma22 = -Conj(s)*delta12 + c*delta22;

                Real* xRealBuf = &XRealBuf[j*ldXReal];
                Real* xImagBuf = &XImagBuf[j*ldXImag];
                // Solve against Q
                const C chi1 = C(xRealBuf[k  ],xImagBuf[k  ]);
                const C chi2 = C(xRealBuf[k+1],xImagBuf[k+1]);
                C eta1 =        c*chi1 + s*chi2;
                C eta2 = -Conj(s)*chi1 + c*chi2;

                // Solve against R
                eta2 /= gamma22;
                eta1 -= gamma12*eta2;
                eta1 /= gamma11;

                // Store the separated real and imaginary comp's of eta{1,2}
                xRealBuf[k+0] = eta1.real();
                xImagBuf[k+0] = eta1.imag();
                xRealBuf[k+1] = eta2.real();
                xImagBuf[k+1] = eta2.imag();

                // Update x0 := x0 - U01 x1
                blas::Axpy
                ( k, -xRealBuf[k  ], &UBuf[ k   *ldu], 1, xRealBuf, 1 );
                blas::Axpy
                ( k, -xImagBuf[k  ], &UBuf[ k   *ldu], 1, xImagBuf, 1 );
                blas::Axpy
                ( k, -xRealBuf[k+1], &UBuf[(k+1)*ldu], 1, xRealBuf, 1 );
                blas::Axpy
                ( k, -xImagBuf[k+1], &UBuf[(k+1)*ldu], 1, xImagBuf, 1 );
            }
        }
        else
        {
            for( Int j=0; j<n; ++j )
            {
                Real* xRealBuf = &XRealBuf[j*ldXReal];
                Real* xImagBuf = &XImagBuf[j*ldXImag];
                // Solve the 1x1 linear system
                C eta1(xRealBuf[k],xImagBuf[k]);
                eta1 /= UBuf[k+k*ldu] - shifts.Get(j,0);
                xRealBuf[k] = eta1.real();
                xImagBuf[k] = eta1.imag();

                // Update x0 := x0 - u01 chi_1
                blas::Axpy( k, -xRealBuf[k], &UBuf[k*ldu], 1, xRealBuf, 1 );
                blas::Axpy( k, -xImagBuf[k], &UBuf[k*ldu], 1, xImagBuf, 1 );
            }
        }
        --k;
    }
}

template<typename F>
inline void
LUN( const Matrix<F>& U, const Matrix<F>& shifts, Matrix<F>& X )
{
    DEBUG_ONLY(CSE cse("msquasitrsm::LUN"))
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

        LUNUnb( U11, shifts, X1 );
        Gemm( NORMAL, NORMAL, F(-1), U01, X1, F(1), X0 );

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
}

template<typename Real>
inline void
LUN
( const Matrix<Real>& U, const Matrix<Complex<Real>>& shifts, 
        Matrix<Real>& XReal, Matrix<Real>& XImag )
{
    DEBUG_ONLY(CSE cse("msquasitrsm::LUN"))
    const Int m = XReal.Height();
    const Int bsize = Blocksize();

    const Int kLast = LastOffset( m, bsize );
    Int k=kLast, kOld=m;
    while( true )
    {
        const bool in2x2 = ( k>0 && U.Get(k,k-1) != Real(0) );
        if( in2x2 )
            --k;
        const Int nb = kOld-k;

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        auto X0Real = XReal( ind0, ALL );
        auto X0Imag = XImag( ind0, ALL );
        auto X1Real = XReal( ind1, ALL );
        auto X1Imag = XImag( ind1, ALL );

        LUNUnb( U11, shifts, X1Real, X1Imag );
        Gemm( NORMAL, NORMAL, Real(-1), U01, X1Real, Real(1), X0Real );
        Gemm( NORMAL, NORMAL, Real(-1), U01, X1Imag, Real(1), X0Imag );

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
}

template<typename F>
inline void
LUNLarge
( const ElementalMatrix<F>& UPre,
  const ElementalMatrix<F>& shiftsPre, 
        ElementalMatrix<F>& XPre )
{
    DEBUG_ONLY(CSE cse("msquasitrsm::LUNLarge"))
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixReadProxy<F,F,VR,STAR> shiftsProx( shiftsPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& X = XProx.Get();
    auto& shifts = shiftsProx.GetLocked();

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
        X1_STAR_VR.AlignWith( shifts );
        X1_STAR_VR    = X1;  // X1[* ,VR] <- X1[MC,MR]
        
        // X1[* ,VR] := U11^-1[* ,* ] X1[* ,VR]
        LocalMultiShiftQuasiTrsm
        ( LEFT, UPPER, NORMAL, F(1), U11_STAR_STAR, shifts, X1_STAR_VR );

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

template<typename Real>
inline void
LUNLarge
( const ElementalMatrix<Real>& UPre, 
  const ElementalMatrix<Complex<Real>>& shiftsPre, 
        ElementalMatrix<Real>& XRealPre, 
        ElementalMatrix<Real>& XImagPre )
{
    DEBUG_ONLY(CSE cse("msquasitrsm::LUNLarge"))
    // TODO: More error checks, especially on alignments?
    typedef Complex<Real> C; 
    const Int m = XRealPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();

    DistMatrixReadProxy<Real,Real,MC,MR> UProx( UPre );
    DistMatrixReadProxy<C,C,VR,STAR> shiftsProx( shiftsPre );
    DistMatrixReadWriteProxy<Real,Real,MC,MR>
      XRealProx( XRealPre ),
      XImagProx( XImagPre );
    auto& U = UProx.GetLocked();
    auto& XReal = XRealProx.Get();
    auto& XImag = XImagProx.Get();
    auto& shifts = shiftsProx.GetLocked();

    DistMatrix<Real,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<Real,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<Real,STAR,MR  > X1Real_STAR_MR(g), X1Imag_STAR_MR(g);
    DistMatrix<Real,STAR,VR  > X1Real_STAR_VR(g), X1Imag_STAR_VR(g);

    const Int kLast = LastOffset( m, bsize );
    Int k=kLast, kOld=m;
    while( true )
    {
        const bool in2x2 = ( k>0 && U.Get(k,k-1) != Real(0) );
        if( in2x2 )
            --k;
        const Int nb = kOld-k;

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        auto X0Real = XReal( ind0, ALL );
        auto X0Imag = XImag( ind0, ALL );
        auto X1Real = XReal( ind1, ALL );
        auto X1Imag = XImag( ind1, ALL );

        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        X1Real_STAR_VR.AlignWith( shifts );
        X1Imag_STAR_VR.AlignWith( shifts );
        X1Real_STAR_VR = X1Real;
        X1Imag_STAR_VR = X1Imag;
        
        // X1[* ,VR] := U11^-1[* ,* ] X1[* ,VR]
        LocalMultiShiftQuasiTrsm
        ( LEFT, UPPER, NORMAL, C(1), 
          U11_STAR_STAR, shifts, X1Real_STAR_VR, X1Imag_STAR_VR );

        X1Real_STAR_MR.AlignWith( X0Real );
        X1Imag_STAR_MR.AlignWith( X0Imag );
        X1Real_STAR_MR = X1Real_STAR_VR; 
        X1Imag_STAR_MR = X1Imag_STAR_VR;
        X1Real = X1Real_STAR_MR;
        X1Imag = X1Imag_STAR_MR;
        U01_MC_STAR.AlignWith( X0Real );
        U01_MC_STAR = U01; 

        // X0[MC,MR] -= U01[MC,* ] X1[* ,MR]
        LocalGemm
        ( NORMAL, NORMAL, 
          Real(-1), U01_MC_STAR, X1Real_STAR_MR, Real(1), X0Real );
        LocalGemm
        ( NORMAL, NORMAL, 
          Real(-1), U01_MC_STAR, X1Imag_STAR_MR, Real(1), X0Imag );

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
}

template<typename F>
inline void
LUNMedium
( const ElementalMatrix<F>& UPre,
  const ElementalMatrix<F>& shiftsPre, 
        ElementalMatrix<F>& XPre )
{
    DEBUG_ONLY(CSE cse("msquasitrsm::LUNMedium"))
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixReadProxy<F,F,VR,STAR> shiftsProx( shiftsPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& X = XProx.Get();
    auto& shifts = shiftsProx.GetLocked();

    DistMatrix<F,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,MR,  STAR> X1Trans_MR_STAR(g);

    DistMatrix<F,MR,  STAR> shifts_MR_STAR(shifts),
                            shifts_MR_STAR_Align(g);

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
        shifts_MR_STAR_Align.AlignWith( X1Trans_MR_STAR );
        shifts_MR_STAR_Align = shifts_MR_STAR;
        LocalMultiShiftQuasiTrsm
        ( RIGHT, UPPER, TRANSPOSE,
          F(1), U11_STAR_STAR, shifts_MR_STAR_Align, X1Trans_MR_STAR );
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

template<typename Real>
inline void
LUNMedium
( const ElementalMatrix<Real>& UPre, 
  const ElementalMatrix<Complex<Real>>& shiftsPre, 
        ElementalMatrix<Real>& XRealPre, 
        ElementalMatrix<Real>& XImagPre )
{
    DEBUG_ONLY(CSE cse("msquasitrsm::LUNMedium"))
    // TODO: Error checks, particularly on alignments?
    typedef Complex<Real> C;
    const Int m = XRealPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();

    DistMatrixReadProxy<Real,Real,MC,MR> UProx( UPre );
    DistMatrixReadProxy<C,C,VR,STAR> shiftsProx( shiftsPre );
    DistMatrixReadWriteProxy<Real,Real,MC,MR>
      XRealProx( XRealPre ),
      XImagProx( XImagPre );
    auto& U = UProx.GetLocked();
    auto& XReal = XRealProx.Get();
    auto& XImag = XImagProx.Get();
    auto& shifts = shiftsProx.GetLocked();

    DistMatrix<Real,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<Real,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<Real,MR,  STAR> X1RealTrans_MR_STAR(g),
                               X1ImagTrans_MR_STAR(g);

    DistMatrix<C,MR,  STAR> shifts_MR_STAR(shifts),
                            shifts_MR_STAR_Align(g);

    const Int kLast = LastOffset( m, bsize );
    Int k=kLast, kOld=m;
    while( true )
    {
        const bool in2x2 = ( k>0 && U.Get(k,k-1) != Real(0) );
        if( in2x2 )
            --k;
        const Int nb = kOld-k;

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        auto X0Real = XReal( ind0, ALL );
        auto X0Imag = XImag( ind0, ALL );
        auto X1Real = XReal( ind1, ALL );
        auto X1Imag = XImag( ind1, ALL );

        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        X1RealTrans_MR_STAR.AlignWith( X0Real );
        X1ImagTrans_MR_STAR.AlignWith( X0Imag );
        Transpose( X1Real, X1RealTrans_MR_STAR );
        Transpose( X1Imag, X1ImagTrans_MR_STAR );
        
        // X1^T[MR,* ] := X1^T[MR,* ] U11^-T[* ,* ]
        //              = (U11^-1[* ,* ] X1[* ,MR])^T
        shifts_MR_STAR_Align.AlignWith( X1RealTrans_MR_STAR );
        shifts_MR_STAR_Align = shifts_MR_STAR;
        LocalMultiShiftQuasiTrsm
        ( RIGHT, UPPER, TRANSPOSE,
          C(1), U11_STAR_STAR, shifts_MR_STAR_Align, 
                X1RealTrans_MR_STAR, X1ImagTrans_MR_STAR );
        Transpose( X1RealTrans_MR_STAR, X1Real );
        Transpose( X1ImagTrans_MR_STAR, X1Imag );

        U01_MC_STAR.AlignWith( X0Real );
        U01_MC_STAR = U01;  // U01[MC,* ] <- U01[MC,MR]

        // X0[MC,MR] -= U01[MC,* ] X1[* ,MR]
        LocalGemm
        ( NORMAL, TRANSPOSE, 
          Real(-1), U01_MC_STAR, X1RealTrans_MR_STAR, Real(1), X0Real );
        LocalGemm
        ( NORMAL, TRANSPOSE, 
          Real(-1), U01_MC_STAR, X1ImagTrans_MR_STAR, Real(1), X0Imag );

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
}

template<typename F,Dist colDist,Dist shiftColDist,Dist shiftRowDist>
inline void
LUNSmall
( const DistMatrix<F,     colDist,STAR        >& U, 
  const DistMatrix<F,shiftColDist,shiftRowDist>& shifts,
        DistMatrix<F,     colDist,STAR        >& X )
{
    DEBUG_ONLY(
      CSE cse("msquasitrsm::LUNSmall");
      AssertSameGrids( U, shifts, X );
      if( U.Height() != U.Width() || U.Width() != X.Height() )
          LogicError
          ("Nonconformal: \n",
           "  U ~ ",U.Height()," x ",U.Width(),"\n",
           "  X ~ ",X.Height()," x ",X.Width(),"\n");
      if( U.ColAlign() != X.ColAlign() )
          LogicError("U and X are assumed to be aligned");
    )
    const Int m = X.Height();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g), X1_STAR_STAR(g),
                            shifts_STAR_STAR(shifts);

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
        LocalMultiShiftQuasiTrsm
        ( LEFT, UPPER, NORMAL, 
          F(1), U11_STAR_STAR, shifts_STAR_STAR, X1_STAR_STAR );
        X1 = X1_STAR_STAR;

        // X0[VC,* ] -= U01[VC,* ] X1[* ,* ]
        LocalGemm( NORMAL, NORMAL, F(-1), U01, X1_STAR_STAR, F(1), X0 );

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
}

template<typename Real,Dist colDist,Dist shiftColDist,Dist shiftRowDist>
inline void
LUNSmall
( const DistMatrix<Real,              colDist,STAR        >& U, 
  const DistMatrix<Complex<Real>,shiftColDist,shiftRowDist>& shifts,
        DistMatrix<Real,              colDist,STAR        >& XReal,
        DistMatrix<Real,              colDist,STAR        >& XImag )
{
    DEBUG_ONLY(
      CSE cse("msquasitrsm::LUNSmall");
      AssertSameGrids( U, shifts, XReal, XImag );
      if( XReal.Height() != XImag.Height() || 
          XReal.Width() != XImag.Width() )
          LogicError("XReal and XImag must be the same size");
      if( U.Height() != U.Width() || U.Width() != XReal.Height() )
          LogicError
          ("Nonconformal: \n",
           "  U ~ ",U.Height()," x ",U.Width(),"\n",
           "  X ~ ",XReal.Height()," x ",XReal.Width(),"\n");
      if( U.ColAlign() != XReal.ColAlign() || 
          U.ColAlign() != XImag.ColAlign() )
          LogicError("U and X are assumed to be aligned");
    )
    typedef Complex<Real> C;
    const Int m = XReal.Height();
    const Int n = XReal.Width();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<Real,STAR,STAR> U11_STAR_STAR(g), X1Real_STAR_STAR(g),
                                                 X1Imag_STAR_STAR(g);
    DistMatrix<C,STAR,STAR> shifts_STAR_STAR(shifts);

    const Int kLast = LastOffset( m, bsize );
    Int k=kLast, kOld=m;
    while( true )
    {
        const bool in2x2 = ( k>0 && U.Get(k,k-1) != Real(0) );
        if( in2x2 )
            --k;
        const Int nb = kOld-k;

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        auto X0Real = XReal( ind0, ALL );
        auto X0Imag = XImag( ind0, ALL );
        auto X1Real = XReal( ind1, ALL );
        auto X1Imag = XImag( ind1, ALL );

        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[VC,* ]
        X1Real_STAR_STAR = X1Real; 
        X1Imag_STAR_STAR = X1Imag; 
        
        // X1[* ,* ] := U11^-1[* ,* ] X1[* ,* ]
        LocalMultiShiftQuasiTrsm
        ( LEFT, UPPER, NORMAL, 
          C(1), U11_STAR_STAR, shifts_STAR_STAR, 
                X1Real_STAR_STAR, X1Imag_STAR_STAR );
        X1Real = X1Real_STAR_STAR;
        X1Imag = X1Imag_STAR_STAR;

        // X0[VC,* ] -= U01[VC,* ] X1[* ,* ]
        LocalGemm
        ( NORMAL, NORMAL, Real(-1), U01, X1Real_STAR_STAR, Real(1), X0Real );
        LocalGemm
        ( NORMAL, NORMAL, Real(-1), U01, X1Imag_STAR_STAR, Real(1), X0Imag );

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
}

} // namespace msquasitrsm
} // namespace El
