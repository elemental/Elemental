/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TRDTRMM_UVAR1_HPP
#define ELEM_TRDTRMM_UVAR1_HPP

#include ELEM_DIAGONALSOLVE_INC
#include ELEM_TRMM_INC

namespace elem {
namespace internal {

template<typename F>
inline void
TrdtrmmUVar1( Matrix<F>& U, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrtdrmmUVar1");
        if( U.Height() != U.Width() )
            LogicError("U must be square");
    )
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    Matrix<F> S01;

    const Int n = U.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto U00 = LockedView( U, 0, 0, k,    k    );
        auto U01 = LockedView( U, 0, k, k,    k+nb );
        auto U11 = LockedView( U, k, k, k+nb, k+nb );
        auto d1 = U11.GetDiagonal();

        S01 = U01;
        DiagonalSolve( LEFT, NORMAL, d1, U01, true );
        Trrk( UPPER, NORMAL, orientation, F(1), U01, S01, F(1), U00 );
        Trmm( RIGHT, UPPER, orientation, UNIT, F(1), U11, U01 );
        TrdtrmmUUnblocked( U11, conjugate );
    }
}

template<typename F>
inline void
TrdtrmmUVar1( DistMatrix<F>& U, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrdtrmmUVar1");
        if( U.Height() != U.Width() )
            LogicError("U must be square");
    )
    const Grid& g = U.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    DistMatrix<F,MC,  STAR> S01_MC_STAR(g);
    DistMatrix<F,VC,  STAR> S01_VC_STAR(g);
    DistMatrix<F,VR,  STAR> U01_VR_STAR(g);
    DistMatrix<F,STAR,MR  > U01Trans_STAR_MR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);

    S01_MC_STAR.AlignWith( U );
    S01_VC_STAR.AlignWith( U );
    U01_VR_STAR.AlignWith( U );
    U01Trans_STAR_MR.AlignWith( U );

    const Int n = U.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto U00 = LockedView( U, 0, 0, k,    k    );
        auto U01 = LockedView( U, 0, k, k,    k+nb );
        auto U11 = LockedView( U, k, k, k+nb, k+nb );
        auto d1 = U11.GetDiagonal();

        S01_MC_STAR = U01;
        S01_VC_STAR = S01_MC_STAR;
        U01_VR_STAR = S01_VC_STAR;
        DiagonalSolve( RIGHT, NORMAL, d1, U01_VR_STAR );
        U01_VR_STAR.TransposePartialColAllGather( U01Trans_STAR_MR, conjugate );
        LocalTrrk( UPPER, F(1), S01_MC_STAR, U01Trans_STAR_MR, F(1), U00 );

        U11_STAR_STAR = U11;
        LocalTrmm
        ( RIGHT, UPPER, orientation, UNIT, F(1), U11_STAR_STAR, U01_VR_STAR );
        U01 = U01_VR_STAR;

        LocalTrdtrmm( UPPER, U11_STAR_STAR, conjugate );
        U11 = U11_STAR_STAR;
    }
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_TRDTRMM_UVAR1_HPP
