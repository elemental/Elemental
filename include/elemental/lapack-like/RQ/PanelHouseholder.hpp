/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_RQ_PANEL_HPP
#define ELEM_LAPACK_RQ_PANEL_HPP

#include "elemental/blas-like/level2/Gemv.hpp"
#include "elemental/blas-like/level2/Ger.hpp"
#include "elemental/lapack-like/Reflector.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {
namespace rq {

template<typename F> 
inline void
PanelHouseholder( Matrix<F>& A, Matrix<F>& t )
{
    DEBUG_ONLY(CallStackEntry cse("rq::PanelHouseholder"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    t.Resize( minDim, 1 );

    const Int iOff = ( n>=m ? 0   : m-n );
    const Int jOff = ( n>=m ? n-m : 0   );

    Matrix<F> z, aBottomRowConj;
    for( Int k=minDim-1; k>=0; --k )
    {
        const Int ki = k + iOff;
        const Int kj = k + jOff;
        auto a10        = ViewRange( A, ki, 0,  ki+1, kj   );
        auto alpha11    = ViewRange( A, ki, kj, ki+1, kj+1 );
        auto ATopPan    = ViewRange( A, 0,  0,  ki,   kj+1 );
        auto aBottomRow = ViewRange( A, ki, 0,  ki+1, kj+1 );

        // Compute the Householder reflector
        const F tau = Reflector( alpha11, a10 );
        t.Set( k, 0, tau );

        // Apply the Householder reflector
        const F alpha = alpha11.Get(0,0);
        alpha11.Set(0,0,1);
        Conjugate( aBottomRow, aBottomRowConj );
        Zeros( z, ATopPan.Height(), 1 );
        Gemv( NORMAL, F(1), ATopPan, aBottomRowConj, F(0), z );
        Ger( -Conj(tau), z, aBottomRowConj, ATopPan );
        alpha11.Set(0,0,alpha);
    }
}

template<typename F> 
inline void
PanelHouseholder( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("rq::PanelHouseholder"))
    Matrix<F> t;
    PanelHouseholder( A, t );
}

template<typename F> 
inline void
PanelHouseholder( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t )
{
    DEBUG_ONLY(
        CallStackEntry cse("rq::PanelHouseholder");
        if( A.Grid() != t.Grid() )
            LogicError("{A,t} must be distributed over the same grid");
        if( !t.AlignedWithDiagonal( A, A.Width()-A.Height() ) )
            LogicError("t must be aligned with A's main diagonal");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    t.Resize( minDim, 1 );

    const Int iOff = ( n>=m ? 0   : m-n );
    const Int jOff = ( n>=m ? n-m : 0   );

    const Grid& g = A.Grid();
    DistMatrix<F> aBottomRowConj(g);
    DistMatrix<F,STAR,MR  > aBottomRowConj_STAR_MR(g);
    DistMatrix<F,MC,  STAR> z_MC_STAR(g);

    for( Int k=minDim-1; k>=0; --k )
    {
        const Int ki = k + iOff;
        const Int kj = k + jOff;
        auto a10        = ViewRange( A, ki, 0,  ki+1, kj   );
        auto alpha11    = ViewRange( A, ki, kj, ki+1, kj+1 );
        auto ATopPan    = ViewRange( A, 0,  0,  ki,   kj+1 );
        auto aBottomRow = ViewRange( A, ki, 0,  ki+1, kj+1 );

        // Compute the Householder reflector
        const F tau = Reflector( alpha11, a10 );
        t.Set( k, 0, tau );

        // Apply the Householder reflector
        const bool myDiagonalEntry = ( g.Row() == alpha11.ColAlign() && 
                                       g.Col() == alpha11.RowAlign() );
        F alpha = 0;
        if( myDiagonalEntry )
        {
            alpha = alpha11.GetLocal(0,0);
            alpha11.SetLocal(0,0,1);
        }
        Conjugate( aBottomRow, aBottomRowConj );
        aBottomRowConj_STAR_MR = aBottomRowConj;
        Zeros( z_MC_STAR, ATopPan.Height(), 1 );
        LocalGemv
        ( NORMAL, F(1), ATopPan, aBottomRowConj_STAR_MR, F(0), z_MC_STAR );
        z_MC_STAR.SumOverRow(); 
        Ger
        ( -Conj(tau), 
          z_MC_STAR.LockedMatrix(),
          aBottomRowConj_STAR_MR.LockedMatrix(),
          ATopPan.Matrix() ); 
        if( myDiagonalEntry )
            alpha11.SetLocal(0,0,alpha);
    }
}

template<typename F> 
inline void
PanelHouseholder( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("rq::PanelHouseholder"))
    DistMatrix<F,MD,STAR> t(A.Grid());
    PanelHouseholder( A, t );
}

} // namespace rq
} // namespace elem

#endif // ifndef ELEM_LAPACK_RQ_PANEL_HPP
