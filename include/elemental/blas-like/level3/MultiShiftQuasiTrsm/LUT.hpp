/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MULTISHIFTQUASITRSM_LUT_HPP
#define ELEM_MULTISHIFTQUASITRSM_LUT_HPP

#include ELEM_GEMM_INC

namespace elem {
namespace msquasitrsm {

template<typename F>
inline void
LUTUnb
( bool conjugate, const Matrix<F>& U, const Matrix<F>& shifts, Matrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("msquasitrsm::LUTUnb"))
    typedef Base<F> Real;
    const Int m = X.Height();
    const Int n = X.Width();
    if( conjugate )
        Conjugate( X );

    const F* UBuf = U.LockedBuffer();
          F* XBuf = X.Buffer();
    const Int ldu = U.LDim();
    const Int ldx = X.LDim();

    Int k=0;
    while( k < m )
    {
        const bool in2x2 = ( k+1<m && UBuf[(k+1)+k*ldu] != F(0) );
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
            const F delta12 = UBuf[ k   +(k+1)*ldu];
            const F delta21 = UBuf[(k+1)+ k   *ldu];
            for( Int j=0; j<n; ++j )
            {
                const F delta11 = UBuf[ k   + k   *ldu] - shifts.Get(j,0);
                const F delta22 = UBuf[(k+1)+(k+1)*ldu] - shifts.Get(j,0);
                // Decompose D = Q R
                Real c; F s;
                const F gamma11 = lapack::Givens( delta11, delta21, &c, &s );
                const F gamma12 =        c*delta12 + s*delta22;
                const F gamma22 = -Conj(s)*delta12 + c*delta22;

                F* xBuf = &XBuf[j*ldx];

                // Solve against R^T
                xBuf[k  ] /= gamma11;
                xBuf[k+1] -= gamma12*xBuf[k];
                xBuf[k+1] /= gamma22;

                // Solve against Q^T
                const F chi1 = xBuf[k  ];
                const F chi2 = xBuf[k+1];
                xBuf[k  ] = c*chi1 - Conj(s)*chi2;
                xBuf[k+1] = s*chi1 +       c*chi2;

                // Update x2 := x2 - U12^T x1
                blas::Axpy
                ( m-(k+2), -xBuf[k  ],
                  &UBuf[ k   +(k+2)*ldu], ldu, &xBuf[k+2], 1 );
                blas::Axpy
                ( m-(k+2), -xBuf[k+1],
                  &UBuf[(k+1)+(k+2)*ldu], ldu, &xBuf[k+2], 1 );
            }
            k += 2;
        }
        else
        {
            for( Int j=0; j<n; ++j )
            {
                F* xBuf = &XBuf[j*ldx];
                xBuf[k] /= UBuf[k+k*ldu] - shifts.Get(j,0);
                blas::Axpy
                ( m-(k+1), -xBuf[k], &UBuf[k+(k+1)*ldu], ldu, &xBuf[k+1], 1 );
            }
            k += 1;
        }
    }
    if( conjugate )
        Conjugate( X );
}

template<typename F>
inline void
LUT
( Orientation orientation, 
  const Matrix<F>& U, const Matrix<F>& shifts, Matrix<F>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("msquasitrsm::LUT");
        if( orientation == NORMAL )
            LogicError("QuasiTrsmLUT expects a (Conjugate)Transpose option");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();

    const bool conjugate = ( orientation==ADJOINT );

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = ( k+nbProp<m && U.Get(k+nbProp,k+nbProp-1) != F(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        auto U11 = LockedViewRange( U, k, k,    k+nb, k+nb );
        auto U12 = LockedViewRange( U, k, k+nb, k+nb, m    );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

        LUTUnb( conjugate, U11, shifts, X1 );
        Gemm( orientation, NORMAL, F(-1), U12, X1, F(1), X2 );
    }
}

// width(X) >> p
template<typename F>
inline void
LUTLarge
( Orientation orientation, 
  const DistMatrix<F>& U, 
  const DistMatrix<F,VR,STAR>& shifts, 
        DistMatrix<F>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("msquasitrsm::LUTLarge");
        if( orientation == NORMAL )
            LogicError("TrsmLUT expects a (Conjugate)Transpose option");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g); 
    DistMatrix<F,STAR,MC  > U12_STAR_MC(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = ( k+nbProp<m && U.Get(k+nbProp,k+nbProp-1) != F(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        auto U11 = LockedViewRange( U, k, k,    k+nb, k+nb );
        auto U12 = LockedViewRange( U, k, k+nb, k+nb, m    );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

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

// width(X) ~= p
template<typename F,Dist shiftColDist,Dist shiftRowDist>
inline void
LUTMedium
( Orientation orientation, 
  const DistMatrix<F>& U, 
  const DistMatrix<F,shiftColDist,shiftRowDist>& shifts, 
        DistMatrix<F>& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("msquasitrsm::LUTMedium");
        if( orientation == NORMAL )
            LogicError("TrsmLUT expects a (Conjugate)Transpose option");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

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

        auto U11 = LockedViewRange( U, k, k,    k+nb, k+nb );
        auto U12 = LockedViewRange( U, k, k+nb, k+nb, m    );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        // X1[* ,VR] <- X1[MC,MR]
        X1Trans_MR_STAR.AlignWith( X2 );
        X1.TransposeColAllGather( X1Trans_MR_STAR, (orientation==ADJOINT) );
        
        // X1[* ,MR] := U11^-[T/H][*,*] X1[* ,MR]
        // X1^[T/H][MR,* ] := X1^[T/H][MR,* ] U11^-1[* ,* ]
        shifts_MR_STAR_Align.AlignWith( X1Trans_MR_STAR );
        shifts_MR_STAR_Align = shifts_MR_STAR;
        LocalMultiShiftQuasiTrsm
        ( RIGHT, UPPER, NORMAL, 
          F(1), U11_STAR_STAR, shifts_MR_STAR_Align, X1Trans_MR_STAR );

        X1.TransposeColFilterFrom( X1Trans_MR_STAR, (orientation==ADJOINT) );
        U12_STAR_MC.AlignWith( X2 );
        U12_STAR_MC = U12; // U12[* ,MC] <- U12[MC,MR]

        // X2[MC,MR] -= (U12[* ,MC])^[T/H] X1[* ,MR]
        //            = U12^[T/H][MC,*] X1[* ,MR]
        LocalGemm
        ( orientation, orientation, 
          F(-1), U12_STAR_MC, X1Trans_MR_STAR, F(1), X2 );
    }
}

// width(X) << p
template<typename F,Dist rowDist,Dist shiftColDist,Dist shiftRowDist>
inline void
LUTSmall
( Orientation orientation, 
  const DistMatrix<F,STAR,             rowDist>& U, 
  const DistMatrix<F,shiftColDist,shiftRowDist>& shifts,
        DistMatrix<F,     rowDist,STAR        >& X )
{
    DEBUG_ONLY(
        CallStackEntry cse("msquasitrsm::LUTSmall");
        if( U.Grid() != X.Grid() )
            LogicError("U and X must be distributed over the same grid");
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
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g), X1_STAR_STAR(g),
                            shifts_STAR_STAR(shifts);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = ( k+nbProp<m && U.Get(k+nbProp,k+nbProp-1) != F(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        auto U11 = LockedViewRange( U, k, k,    k+nb, k+nb );
        auto U12 = LockedViewRange( U, k, k+nb, k+nb, m    );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

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

} // namespace msquasitrsm
} // namespace elem

#endif // ifndef ELEM_MULTISHIFTQUASITRSM_LUT_HPP
