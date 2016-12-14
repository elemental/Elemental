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
void LUTUnb
( bool conjugate, const Matrix<F>& U, const Matrix<F>& shifts, Matrix<F>& X )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Int m = X.Height();
    const Int n = X.Width();

    if( conjugate )
        Conjugate( X );

    const F* UBuf = U.LockedBuffer();
          F* XBuf = X.Buffer();
    const Int ldU = U.LDim();
    const Int ldX = X.LDim();

    Int k=0;
    while( k < m )
    {
        const bool in2x2 = ( k+1<m && UBuf[(k+1)+k*ldU] != F(0) );
        if( in2x2 )
        {
            // Solve the 2x2 linear systems via 2x2 QR decompositions produced
            // by the Givens rotation
            //    | c        s | | U(k,  k)-shift | = | gamma11 | 
            //    | -conj(s) c | | U(k+1,k)       |   | 0       |
            //
            // and by also forming the right two entries of the 2x2 resulting
            // upper-triangular matrix, say gamma12 and gamma22
            //
            // Extract the constant part of the 2x2 diagonal block, D
            const F delta12 = UBuf[ k   +(k+1)*ldU];
            const F delta21 = UBuf[(k+1)+ k   *ldU];
            for( Int j=0; j<n; ++j )
            {
                const F delta11 = UBuf[ k   + k   *ldU] - shifts.Get(j,0);
                const F delta22 = UBuf[(k+1)+(k+1)*ldU] - shifts.Get(j,0);
                // Decompose D = Q R
                Real c; F s;
                const F gamma11 = Givens( delta11, delta21, c, s );
                const F gamma12 =        c*delta12 + s*delta22;
                const F gamma22 = -Conj(s)*delta12 + c*delta22;

                F* xBuf = &XBuf[j*ldX];

                // Solve against R^T
                F chi1 = xBuf[k  ];
                F chi2 = xBuf[k+1];
                chi1 /= gamma11;
                chi2 -= gamma12*chi1;
                chi2 /= gamma22;

                // Solve against Q^T
                xBuf[k  ] = c*chi1 - Conj(s)*chi2;
                xBuf[k+1] = s*chi1 +       c*chi2;

                // Update x2 := x2 - U12^T x1
                blas::Axpy
                ( m-(k+2), -xBuf[k  ],
                  &UBuf[ k   +(k+2)*ldU], ldU, &xBuf[k+2], 1 );
                blas::Axpy
                ( m-(k+2), -xBuf[k+1],
                  &UBuf[(k+1)+(k+2)*ldU], ldU, &xBuf[k+2], 1 );
            }
            k += 2;
        }
        else
        {
            for( Int j=0; j<n; ++j )
            {
                F* xBuf = &XBuf[j*ldX];
                xBuf[k] /= UBuf[k+k*ldU] - shifts.Get(j,0);
                blas::Axpy
                ( m-(k+1), -xBuf[k], &UBuf[k+(k+1)*ldU], ldU, &xBuf[k+1], 1 );
            }
            k += 1;
        }
    }
    if( conjugate )
        Conjugate( X );
}

template<typename Real>
void LUTUnb
( bool conjugate, 
  const Matrix<Real>& U, 
  const Matrix<Complex<Real>>& shifts, 
        Matrix<Real>& XReal,
        Matrix<Real>& XImag )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;
    const Int m = XReal.Height();
    const Int n = XReal.Width();
  
    if( conjugate )
        XImag *= -1;

    const Real* UBuf = U.LockedBuffer();
          Real* XRealBuf = XReal.Buffer();
          Real* XImagBuf = XImag.Buffer();
    const Int ldU = U.LDim();
    const Int ldXReal = XReal.LDim();
    const Int ldXImag = XImag.LDim();

    Int k=0;
    while( k < m )
    {
        const bool in2x2 = ( k+1<m && UBuf[(k+1)+k*ldU] != Real(0) );
        if( in2x2 )
        {
            // Solve the 2x2 linear systems via 2x2 QR decompositions produced
            // by the Givens rotation
            //    | c        s | | U(k,  k)-shift | = | gamma11 | 
            //    | -conj(s) c | | U(k+1,k)       |   | 0       |
            //
            // and by also forming the right two entries of the 2x2 resulting
            // upper-triangular matrix, say gamma12 and gamma22
            //
            // Extract the constant part of the 2x2 diagonal block, D
            const Real delta12 = UBuf[ k   +(k+1)*ldU];
            const Real delta21 = UBuf[(k+1)+ k   *ldU];
            for( Int j=0; j<n; ++j )
            {
                const C delta11 = UBuf[ k   + k   *ldU] - shifts.Get(j,0);
                const C delta22 = UBuf[(k+1)+(k+1)*ldU] - shifts.Get(j,0);
                // Decompose D = Q R
                Real c; C s;
                const C gamma11 = Givens( delta11, C(delta21), c, s );
                const C gamma12 =        c*delta12 + s*delta22;
                const C gamma22 = -Conj(s)*delta12 + c*delta22;

                Real* xRealBuf = &XRealBuf[j*ldXReal];
                Real* xImagBuf = &XImagBuf[j*ldXImag]; 

                // Solve against R^T
                C chi1(xRealBuf[k  ],xImagBuf[k  ]);
                C chi2(xRealBuf[k+1],xImagBuf[k+1]);
                chi1 /= gamma11;
                chi2 -= gamma12*chi1;
                chi2 /= gamma22;

                // Solve against Q^T
                const C eta1 = c*chi1 - Conj(s)*chi2;
                const C eta2 = s*chi1 +       c*chi2;
                xRealBuf[k  ] = eta1.real();
                xImagBuf[k  ] = eta1.imag();
                xRealBuf[k+1] = eta2.real();
                xImagBuf[k+1] = eta2.imag();

                // Update x2 := x2 - U12^T x1
                blas::Axpy
                ( m-(k+2), -xRealBuf[k  ],
                  &UBuf[ k   +(k+2)*ldU], ldU, &xRealBuf[k+2], 1 );
                blas::Axpy
                ( m-(k+2), -xImagBuf[k  ],
                  &UBuf[ k   +(k+2)*ldU], ldU, &xImagBuf[k+2], 1 );
                blas::Axpy
                ( m-(k+2), -xRealBuf[k+1],
                  &UBuf[(k+1)+(k+2)*ldU], ldU, &xRealBuf[k+2], 1 );
                blas::Axpy
                ( m-(k+2), -xImagBuf[k+1],
                  &UBuf[(k+1)+(k+2)*ldU], ldU, &xImagBuf[k+2], 1 );
            }
            k += 2;
        }
        else
        {
            for( Int j=0; j<n; ++j )
            {
                Real* xRealBuf = &XRealBuf[j*ldXReal];
                Real* xImagBuf = &XImagBuf[j*ldXImag];
                C eta1( xRealBuf[k], xImagBuf[k] );
                eta1 /= UBuf[k+k*ldU] - shifts.Get(j,0);
                xRealBuf[k] = eta1.real();
                xImagBuf[k] = eta1.imag();
                blas::Axpy
                ( m-(k+1), -xRealBuf[k], 
                  &UBuf[k+(k+1)*ldU], ldU, &xRealBuf[k+1], 1 );
                blas::Axpy
                ( m-(k+1), -xImagBuf[k], 
                  &UBuf[k+(k+1)*ldU], ldU, &xImagBuf[k+1], 1 );
            }
            k += 1;
        }
    }
    if( conjugate )
        XImag *= -1;
}

template<typename F>
void LUT
( Orientation orientation, 
  const Matrix<F>& U, const Matrix<F>& shifts, Matrix<F>& X )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( orientation == NORMAL )
          LogicError("QuasiTrsmLUT expects a (Conjugate)Transpose option");
    )
    const Int m = X.Height();
    const Int bsize = Blocksize();

    const bool conjugate = ( orientation==ADJOINT );
    if( conjugate )
        Conjugate( X );

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = ( k+nbProp<m && U.Get(k+nbProp,k+nbProp-1) != F(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, m    );

        auto U11 = U( ind1, ind1 );
        auto U12 = U( ind1, ind2 );

        auto X1 = X( ind1, ALL );
        auto X2 = X( ind2, ALL );

        LUTUnb( false, U11, shifts, X1 );
        Gemm( TRANSPOSE, NORMAL, F(-1), U12, X1, F(1), X2 );
    }

    if( conjugate )
        Conjugate( X );
}

template<typename Real>
void LUT
( Orientation orientation, 
  const Matrix<Real>& U, 
  const Matrix<Complex<Real>>& shifts, 
        Matrix<Real>& XReal, Matrix<Real>& XImag )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( orientation == NORMAL )
          LogicError("QuasiTrsmLUT expects a (Conjugate)Transpose option");
    )
    const Int m = XReal.Height();
    const Int bsize = Blocksize();

    const bool conjugate = ( orientation==ADJOINT );
    if( conjugate )
        XImag *= -1;

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = 
            ( k+nbProp<m && U.Get(k+nbProp,k+nbProp-1) != Real(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, m    );

        auto U11 = U( ind1, ind1 );
        auto U12 = U( ind1, ind2 );

        auto X1Real = XReal( ind1, ALL );
        auto X1Imag = XImag( ind1, ALL );
        auto X2Real = XReal( ind2, ALL );
        auto X2Imag = XImag( ind2, ALL );

        LUTUnb( false, U11, shifts, X1Real, X1Imag );
        Gemm( TRANSPOSE, NORMAL, Real(-1), U12, X1Real, Real(1), X2Real );
        Gemm( TRANSPOSE, NORMAL, Real(-1), U12, X1Imag, Real(1), X2Imag );
    }
    if( conjugate )
        XImag *= -1;
}

// width(X) >> p
template<typename F>
void LUTLarge
( Orientation orientation, 
  const AbstractDistMatrix<F>& UPre,
  const AbstractDistMatrix<F>& shiftsPre,
        AbstractDistMatrix<F>& XPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( orientation == NORMAL )
          LogicError("TrsmLUT expects a (Conjugate)Transpose option");
    )
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixReadProxy<F,F,VR,STAR> shiftsProx( shiftsPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& X = XProx.Get();
    auto& shifts = shiftsProx.GetLocked();

    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g); 
    DistMatrix<F,STAR,MC  > U12_STAR_MC(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = ( k+nbProp<m && U.Get(k+nbProp,k+nbProp-1) != F(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, m    );

        auto U11 = U( ind1, ind1 );
        auto U12 = U( ind1, ind2 );

        auto X1 = X( ind1, ALL );
        auto X2 = X( ind2, ALL );

        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        X1_STAR_VR.AlignWith( shifts );
        X1_STAR_VR    = X1;  // X1[* ,VR] <- X1[MC,MR]
        
        // X1[* ,VR] := U11^-[T/H][*,*] X1[* ,VR]
        LocalMultiShiftQuasiTrsm
        ( LEFT, UPPER, orientation, F(1), U11_STAR_STAR, shifts, X1_STAR_VR );

        X1_STAR_MR.AlignWith( X2 );
        X1_STAR_MR  = X1_STAR_VR; // X1[* ,MR]  <- X1[* ,VR]
        X1          = X1_STAR_MR; // X1[MC,MR]  <- X1[* ,MR]
        U12_STAR_MC.AlignWith( X2 );
        U12_STAR_MC = U12;        // U12[* ,MC] <- U12[MC,MR]

        // X2[MC,MR] -= (U12[* ,MC])^(T/H) X1[* ,MR]
        //            = U12^(T/H)[MC,*] X1[* ,MR]
        LocalGemm
        ( orientation, NORMAL, F(-1), U12_STAR_MC, X1_STAR_MR, F(1), X2 );
    }
}

template<typename Real>
void LUTLarge
( Orientation orientation, 
  const AbstractDistMatrix<Real>& UPre, 
  const AbstractDistMatrix<Complex<Real>>& shiftsPre, 
        AbstractDistMatrix<Real>& XRealPre, 
        AbstractDistMatrix<Real>& XImagPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( orientation == NORMAL )
          LogicError("TrsmLUT expects a (Conjugate)Transpose option");
    )
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

    DistMatrix<Real,STAR,STAR> U11_STAR_STAR(g); 
    DistMatrix<Real,STAR,MC  > U12_STAR_MC(g);
    DistMatrix<Real,STAR,MR  > X1Real_STAR_MR(g), X1Imag_STAR_MR(g);
    DistMatrix<Real,STAR,VR  > X1Real_STAR_VR(g), X1Imag_STAR_VR(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = 
            ( k+nbProp<m && U.Get(k+nbProp,k+nbProp-1) != Real(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, m    );

        auto U11 = U( ind1, ind1 );
        auto U12 = U( ind1, ind2 );

        auto X1Real = XReal( ind1, ALL );
        auto X1Imag = XImag( ind1, ALL );
        auto X2Real = XReal( ind2, ALL );
        auto X2Imag = XImag( ind2, ALL );

        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        X1Real_STAR_VR.AlignWith( shifts );
        X1Imag_STAR_VR.AlignWith( shifts );
        X1Real_STAR_VR = X1Real; 
        X1Imag_STAR_VR = X1Imag; 
        
        // X1[* ,VR] := U11^-[T/H][*,*] X1[* ,VR]
        LocalMultiShiftQuasiTrsm
        ( LEFT, UPPER, orientation, 
          C(1), U11_STAR_STAR, shifts, X1Real_STAR_VR, X1Imag_STAR_VR );

        X1Real_STAR_MR.AlignWith( X2Real );
        X1Imag_STAR_MR.AlignWith( X2Imag );
        X1Real_STAR_MR = X1Real_STAR_VR; 
        X1Imag_STAR_MR = X1Imag_STAR_VR; 
        X1Real = X1Real_STAR_MR;
        X1Imag = X1Imag_STAR_MR;
        U12_STAR_MC.AlignWith( X2Real );
        U12_STAR_MC = U12; 

        // X2[MC,MR] -= (U12[* ,MC])^(T/H) X1[* ,MR]
        //            = U12^(T/H)[MC,*] X1[* ,MR]
        LocalGemm
        ( orientation, NORMAL, 
          Real(-1), U12_STAR_MC, X1Real_STAR_MR, Real(1), X2Real );
        LocalGemm
        ( orientation, NORMAL, 
          Real(-1), U12_STAR_MC, X1Imag_STAR_MR, Real(1), X2Imag );
    }
}

// width(X) ~= p
template<typename F>
void LUTMedium
( Orientation orientation, 
  const AbstractDistMatrix<F>& UPre,
  const AbstractDistMatrix<F>& shiftsPre, 
        AbstractDistMatrix<F>& XPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( orientation == NORMAL )
          LogicError("TrsmLUT expects a (Conjugate)Transpose option");
    )
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixReadProxy<F,F,VR,STAR> shiftsProx( shiftsPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& X = XProx.Get();
    auto& shifts = shiftsProx.GetLocked();

    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g); 
    DistMatrix<F,STAR,MC  > U12_STAR_MC(g);
    DistMatrix<F,MR,  STAR> X1Trans_MR_STAR(g);

    DistMatrix<F,MR,  STAR> shifts_MR_STAR(shifts),
                            shifts_MR_STAR_Align(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = ( k+nbProp<m && U.Get(k+nbProp,k+nbProp-1) != F(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, m    );

        auto U11 = U( ind1, ind1 );
        auto U12 = U( ind1, ind2 );

        auto X1 = X( ind1, ALL );
        auto X2 = X( ind2, ALL );

        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        // X1[* ,VR] <- X1[MC,MR]
        X1Trans_MR_STAR.AlignWith( X2 );
        Transpose( X1, X1Trans_MR_STAR, (orientation==ADJOINT) );
        
        // X1[* ,MR] := U11^-[T/H][*,*] X1[* ,MR]
        // X1^[T/H][MR,* ] := X1^[T/H][MR,* ] U11^-1[* ,* ]
        shifts_MR_STAR_Align.AlignWith( X1Trans_MR_STAR );
        shifts_MR_STAR_Align = shifts_MR_STAR;
        LocalMultiShiftQuasiTrsm
        ( RIGHT, UPPER, NORMAL, 
          F(1), U11_STAR_STAR, shifts_MR_STAR_Align, X1Trans_MR_STAR );

        Transpose( X1Trans_MR_STAR, X1, (orientation==ADJOINT) );
        U12_STAR_MC.AlignWith( X2 );
        U12_STAR_MC = U12; // U12[* ,MC] <- U12[MC,MR]

        // X2[MC,MR] -= (U12[* ,MC])^[T/H] X1[* ,MR]
        //            = U12^[T/H][MC,*] X1[* ,MR]
        LocalGemm
        ( orientation, orientation, 
          F(-1), U12_STAR_MC, X1Trans_MR_STAR, F(1), X2 );
    }
}

template<typename Real>
void LUTMedium
( Orientation orientation, 
  const AbstractDistMatrix<Real>& UPre, 
  const AbstractDistMatrix<Complex<Real>>& shiftsPre, 
        AbstractDistMatrix<Real>& XRealPre, 
        AbstractDistMatrix<Real>& XImagPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( orientation == NORMAL )
          LogicError("TrsmLUT expects a (Conjugate)Transpose option");
    )
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

    DistMatrix<Real,STAR,STAR> U11_STAR_STAR(g); 
    DistMatrix<Real,STAR,MC  > U12_STAR_MC(g);
    DistMatrix<Real,MR,  STAR> X1RealTrans_MR_STAR(g), X1ImagTrans_MR_STAR(g);

    DistMatrix<Complex<Real>,MR,STAR> shifts_MR_STAR(shifts),
                                      shifts_MR_STAR_Align(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = 
            ( k+nbProp<m && U.Get(k+nbProp,k+nbProp-1) != Real(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, m    );

        auto U11 = U( ind1, ind1 );
        auto U12 = U( ind1, ind2 );

        auto X1Real = XReal( ind1, ALL );
        auto X1Imag = XImag( ind1, ALL );
        auto X2Real = XReal( ind2, ALL );
        auto X2Imag = XImag( ind2, ALL );

        U11_STAR_STAR = U11; 
        X1RealTrans_MR_STAR.AlignWith( X2Real );
        X1ImagTrans_MR_STAR.AlignWith( X2Imag );
        Transpose( X1Real, X1RealTrans_MR_STAR, (orientation==ADJOINT) );
        Transpose( X1Imag, X1ImagTrans_MR_STAR, (orientation==ADJOINT) );
        
        // X1[* ,MR] := U11^-[T/H][*,*] X1[* ,MR]
        // X1^[T/H][MR,* ] := X1^[T/H][MR,* ] U11^-1[* ,* ]
        shifts_MR_STAR_Align.AlignWith( X1RealTrans_MR_STAR );
        shifts_MR_STAR_Align = shifts_MR_STAR;
        LocalMultiShiftQuasiTrsm
        ( RIGHT, UPPER, NORMAL, 
          C(1), U11_STAR_STAR, shifts_MR_STAR_Align, 
                X1RealTrans_MR_STAR, X1ImagTrans_MR_STAR );

        Transpose( X1RealTrans_MR_STAR, X1Real, (orientation==ADJOINT) );
        Transpose( X1ImagTrans_MR_STAR, X1Imag, (orientation==ADJOINT) );
        U12_STAR_MC.AlignWith( X2Real );
        U12_STAR_MC = U12;

        // X2[MC,MR] -= (U12[* ,MC])^[T/H] X1[* ,MR]
        //            = U12^[T/H][MC,*] X1[* ,MR]
        LocalGemm
        ( orientation, orientation, 
          Real(-1), U12_STAR_MC, X1RealTrans_MR_STAR, Real(1), X2Real );
        LocalGemm
        ( orientation, orientation, 
          Real(-1), U12_STAR_MC, X1ImagTrans_MR_STAR, Real(1), X2Imag );
    }
}

// width(X) << p
template<typename F,Dist rowDist,Dist shiftColDist,Dist shiftRowDist>
void LUTSmall
( Orientation orientation, 
  const DistMatrix<F,STAR,             rowDist>& U, 
  const DistMatrix<F,shiftColDist,shiftRowDist>& shifts,
        DistMatrix<F,     rowDist,STAR        >& X )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( U, shifts, X );
      if( orientation == NORMAL )
          LogicError("TrsmLUT expects a (Conjugate)Transpose option");
      if( U.Height() != U.Width() || U.Height() != X.Height() )
          LogicError
          ("Nonconformal: \n",
           "  U ~ ",U.Height()," x ",U.Width(),"\n",
           "  X ~ ",X.Height()," x ",X.Width(),"\n");
      if( U.RowAlign() != X.ColAlign() )
          LogicError("U and X are assumed to be aligned");
    )
    const Int m = X.Height();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g), X1_STAR_STAR(g),
                            shifts_STAR_STAR(shifts);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = ( k+nbProp<m && U.Get(k+nbProp,k+nbProp-1) != F(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, m    );

        auto U11 = U( ind1, ind1 );
        auto U12 = U( ind1, ind2 );

        auto X1 = X( ind1, ALL );
        auto X2 = X( ind2, ALL );

        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[* ,VR]
        X1_STAR_STAR = X1;   // X1[* ,* ] <- X1[VR,* ]
        
        // X1[* ,* ] := U11^-[T/H][* ,* ] X1[* ,* ]
        LocalMultiShiftQuasiTrsm
        ( LEFT, UPPER, orientation,
          F(1), U11_STAR_STAR, shifts_STAR_STAR, X1_STAR_STAR );

        X1 = X1_STAR_STAR;

        // X2[VR,* ] -= U12[* ,VR]^[T/H] X1[* ,* ]
        LocalGemm( orientation, NORMAL, F(-1), U12, X1_STAR_STAR, F(1), X2 );
    }
}

template<typename Real,Dist rowDist,Dist shiftColDist,Dist shiftRowDist>
void LUTSmall
( Orientation orientation, 
  const DistMatrix<Real,STAR,         rowDist>& U, 
  const DistMatrix<Complex<Real>,shiftColDist,shiftRowDist>& shifts,
        DistMatrix<Real,              rowDist,STAR        >& XReal,
        DistMatrix<Real,              rowDist,STAR        >& XImag )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( U, shifts, XReal, XImag );
      if( orientation == NORMAL )
          LogicError("TrsmLUT expects a (Conjugate)Transpose option");
      if( U.Height() != U.Width() || U.Height() != XReal.Height() )
          LogicError
          ("Nonconformal: \n",
           "  U ~ ",U.Height()," x ",U.Width(),"\n",
           "  X ~ ",XReal.Height()," x ",XReal.Width(),"\n");
      if( U.RowAlign() != XReal.ColAlign() ||
          U.RowAlign() != XImag.ColAlign() )
          LogicError("U and X are assumed to be aligned");
    )
    typedef Complex<Real> C;
    const Int m = XReal.Height();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<Real,STAR,STAR> U11_STAR_STAR(g), X1Real_STAR_STAR(g),
                                                 X1Imag_STAR_STAR(g);
    DistMatrix<C,STAR,STAR> shifts_STAR_STAR(shifts);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = 
            ( k+nbProp<m && U.Get(k+nbProp,k+nbProp-1) != Real(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, m    );

        auto U11 = U( ind1, ind1 );
        auto U12 = U( ind1, ind2 );

        auto X1Real = XReal( ind1, ALL );
        auto X1Imag = XImag( ind1, ALL );
        auto X2Real = XReal( ind2, ALL );
        auto X2Imag = XImag( ind2, ALL );

        U11_STAR_STAR = U11;
        X1Real_STAR_STAR = X1Real;  
        X1Imag_STAR_STAR = X1Imag;  
        
        // X1[* ,* ] := U11^-[T/H][* ,* ] X1[* ,* ]
        LocalMultiShiftQuasiTrsm
        ( LEFT, UPPER, orientation,
          C(1), U11_STAR_STAR, shifts_STAR_STAR, 
                X1Real_STAR_STAR, X1Imag_STAR_STAR );

        X1Real = X1Real_STAR_STAR;
        X1Imag = X1Imag_STAR_STAR;

        // X2[VR,* ] -= U12[* ,VR]^[T/H] X1[* ,* ]
        LocalGemm
        ( orientation, NORMAL, 
          Real(-1), U12, X1Real_STAR_STAR, Real(1), X2Real );
        LocalGemm
        ( orientation, NORMAL, 
          Real(-1), U12, X1Imag_STAR_STAR, Real(1), X2Imag );
    }
}

} // namespace msquasitrsm
} // namespace El
