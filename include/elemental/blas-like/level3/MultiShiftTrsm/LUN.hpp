/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MULTISHIFTTRSM_LUN_HPP
#define ELEM_MULTISHIFTTRSM_LUN_HPP

namespace elem {
namespace mstrsm {

template<typename F>
inline void
LeftUnb
( UpperOrLower uplo, Orientation orientation, F alpha, 
  Matrix<F>& T, const Matrix<F>& shifts, Matrix<F>& X ) 
{
    DEBUG_ONLY(
        CallStackEntry cse("mstrsm::LeftUnb");
        if( shifts.Height() != X.Width() )
            LogicError("Incompatible number of shifts");
    )
    const char uploChar = ( uplo==LOWER ? 'L' : 'U' );
    char orientChar; 
    switch( orientation )
    {
    case NORMAL:    orientChar = 'N'; break;
    case TRANSPOSE: orientChar = 'T'; break;
    default:        orientChar = 'C';
    }
    auto diag = T.GetDiagonal();
    const Int n = T.Height();
    const Int ldim = T.LDim();
    const Int numShifts = shifts.Height();
    Scale( alpha, X );
    for( Int j=0; j<numShifts; ++j )
    {
        UpdateDiagonal( T, -shifts.Get(j,0) );
        blas::Trsv
        ( uploChar, orientChar, 'N', n, 
          T.LockedBuffer(), ldim, X.Buffer(0,j), 1 );
        T.SetDiagonal( diag );
    }
}

template<typename F>
inline void
LUN( F alpha, Matrix<F>& U, const Matrix<F>& shifts, Matrix<F>& X ) 
{
    DEBUG_ONLY(CallStackEntry cse("mstrsm::LUN"))
    Scale( alpha, X );
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( m, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto U01 = LockedViewRange( U, 0, k, k,    k+nb );
        auto U11 =       ViewRange( U, k, k, k+nb, k+nb );

        auto X0 = ViewRange( X, 0, 0, k,    n );
        auto X1 = ViewRange( X, k, 0, k+nb, n );

        LeftUnb( UPPER, NORMAL, F(1), U11, shifts, X1 );
        Gemm( NORMAL, NORMAL, F(-1), U01, X1, F(1), X0 );
    }
}

template<typename F>
inline void
LUN
( F alpha, const DistMatrix<F>& U, const DistMatrix<F,VR,STAR>& shifts,
  DistMatrix<F>& X ) 
{
    DEBUG_ONLY(CallStackEntry cse("mstrsm::LUN"))
    Scale( alpha, X );

    const Grid& g = U.Grid();
    DistMatrix<F,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( m, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto U01 = LockedViewRange( U, 0, k, k,    k+nb );
        auto U11 = LockedViewRange( U, k, k, k+nb, k+nb );

        auto X0 = ViewRange( X, 0, 0, k,    n );
        auto X1 = ViewRange( X, k, 0, k+nb, n );

        // X1[* ,VR] := U11^-1[* ,* ] X1[* ,VR]
        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        X1_STAR_VR.AlignWith( shifts );
        X1_STAR_VR = X1; // X1[* ,VR] <- X1[MC,MR]
        LUN
        ( F(1), U11_STAR_STAR.Matrix(), shifts.LockedMatrix(), 
          X1_STAR_VR.Matrix() );

        X1_STAR_MR.AlignWith( X0 );
        X1_STAR_MR = X1_STAR_VR; // X1[* ,MR]  <- X1[* ,VR]
        X1 = X1_STAR_MR; // X1[MC,MR] <- X1[* ,MR]

        // X0[MC,MR] -= U01[MC,* ] X1[* ,MR]
        U01_MC_STAR.AlignWith( X0 );
        U01_MC_STAR = U01; // U01[MC,* ] <- U01[MC,MR]
        LocalGemm( NORMAL, NORMAL, F(-1), U01_MC_STAR, X1_STAR_MR, F(1), X0 );
    }
}

} // namespace mstrsm
} // namespace elem

#endif // ifndef ELEM_MULTISHIFTTRSM_LUN_HPP
