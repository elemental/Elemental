/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_LQ_PANELHOUSEHOLDER_HPP
#define ELEM_LAPACK_LQ_PANELHOUSEHOLDER_HPP

#include "elemental/blas-like/level1/Conjugate.hpp"
#include "elemental/blas-like/level2/Gemv.hpp"
#include "elemental/blas-like/level2/Ger.hpp"
#include "elemental/lapack-like/Reflector.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {
namespace lq {

template<typename F>
inline void
PanelHouseholder( Matrix<F>& A, Matrix<F>& t )
{
    DEBUG_ONLY(CallStackEntry cse("lq::PanelHouseholder"))
    Matrix<F> z, aTopRowConj;

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    DEBUG_ONLY(
        if( t.Height() != minDim || t.Width() != 1 )
            LogicError("t must be a vector of length minDim(A)");
    )
    for( Int k=0; k<minDim; ++k )
    {
        auto alpha11    = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12        = ViewRange( A, k,   k+1, k+1, n   );
        auto aTopRow    = ViewRange( A, k,   k,   k+1, n   );
        auto ABottomPan = ViewRange( A, k+1, k,   m,   n   );

        // Compute the Householder reflector
        const F tau = Reflector( alpha11, a12 );
        t.Set( k, 0, tau );

        // Apply the Householder reflector
        const F alpha = alpha11.Get(0,0);
        alpha11.Set(0,0,1);
        Conjugate( aTopRow, aTopRowConj );
        Zeros( z, ABottomPan.Height(), 1 );
        Gemv( NORMAL, F(1), ABottomPan, aTopRowConj, F(0), z );
        Ger( -Conj(tau), z, aTopRowConj, ABottomPan );
        alpha11.Set(0,0,alpha);
    }
}

template<typename F>
inline void
PanelHouseholder( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("lq::PanelHouseholder"))
    Matrix<F> t;
    PanelHouseholder( A, t );
}

template<typename F>
inline void
PanelHouseholder( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t )
{
    DEBUG_ONLY(
        CallStackEntry cse("lq::PanelHouseholder");
        if( A.Grid() != t.Grid() )
            LogicError("{A,t} must be distributed over the same grid");
        if( !t.AlignedWithDiagonal( A, 0 ) )
            LogicError("t must be aligned with A's main diagonal");
    )
    const Grid& g = A.Grid();
    DistMatrix<F> aTopRowConj(g);
    DistMatrix<F,STAR,MR  > aTopRowConj_STAR_MR(g);
    DistMatrix<F,MC,  STAR> z_MC_STAR(g);

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    DEBUG_ONLY(
        if( t.Height() != minDim || t.Width() != 1 )
            LogicError("t must be a vector of length minDim(A)");
    )
    for( Int k=0; k<minDim; ++k )
    {
        auto alpha11    = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12        = ViewRange( A, k,   k+1, k+1, n   );
        auto aTopRow    = ViewRange( A, k,   k,   k+1, n   );
        auto ABottomPan = ViewRange( A, k+1, k,   m,   n   );

        // Compute the Householder reflector
        const F tau = Reflector( alpha11, a12 );
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
        aTopRowConj_STAR_MR.AlignWith( ABottomPan );
        Conjugate( aTopRow, aTopRowConj );
        aTopRowConj_STAR_MR = aTopRowConj;
        z_MC_STAR.AlignWith( ABottomPan );
        Zeros( z_MC_STAR, ABottomPan.Height(), 1 );
        LocalGemv
        ( NORMAL, F(1), ABottomPan, aTopRowConj_STAR_MR, F(0), z_MC_STAR );
        z_MC_STAR.SumOverRow();
        Ger
        ( -Conj(tau),
          z_MC_STAR.LockedMatrix(),
          aTopRowConj_STAR_MR.LockedMatrix(),
          ABottomPan.Matrix() );
        if( myDiagonalEntry )
            alpha11.SetLocal(0,0,alpha);
    }
}

template<typename F>
inline void
PanelHouseholder( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("lq::PanelHouseholder"))
    DistMatrix<F,MD,STAR> t(A.Grid());
    PanelHouseholder( A, t );
}

} // namespace lq
} // namespace elem

#endif // ifndef ELEM_LAPACK_LQ_PANELHOUSEHOLDER_HPP
