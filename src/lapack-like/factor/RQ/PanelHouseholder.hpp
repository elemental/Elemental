/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_RQ_PANEL_HPP
#define EL_RQ_PANEL_HPP



namespace El {
namespace rq {

template<typename F> 
inline void
PanelHouseholder( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d )
{
    DEBUG_ONLY(CallStackEntry cse("rq::PanelHouseholder"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int iOff = m-minDim;
    const Int jOff = n-minDim;

    t.Resize( minDim, 1 );
    d.Resize( minDim, 1 );

    Matrix<F> z01;
    for( Int k=minDim-1; k>=0; --k )
    {
        const Int ki = k + iOff;
        const Int kj = k + jOff;
        auto a10     = ViewRange( A, ki, 0,  ki+1, kj   );
        auto alpha11 = ViewRange( A, ki, kj, ki+1, kj+1 );
        auto A0L     = ViewRange( A, 0,  0,  ki,   kj+1 );
        auto a1L     = ViewRange( A, ki, 0,  ki+1, kj+1 );

        // Find tau and v such that
        //  |a10 alpha11| /I - tau |v^T| |conj(v) 1|\ = |0 beta|
        //                \        |1  |            /
        const F tau = RightReflector( alpha11, a10 );
        t.Set( k, 0, tau );

        // Temporarily set a1L = | v 1 |
        const F alpha = alpha11.Get(0,0);
        alpha11.Set(0,0,1);

        // A2R := A2R Hous(a1L^T,tau)
        //      = A2R (I - tau a1L^T conj(a1L))
        //      = A2R - tau (A2R a1L^T) conj(a1L)
        Zeros( z01, A0L.Height(), 1 );
        Gemv( NORMAL, F(1), A0L, a1L, F(0), z01 );
        Ger( -tau, z01, a1L, A0L );

        // Reset alpha11's value
        alpha11.Set(0,0,alpha);
    }
    // Form d and rescale R
    auto R = View( A, 0, jOff, m, minDim );
    d = R.GetRealPartOfDiagonal();
    typedef Base<F> Real;
    for( Int j=0; j<minDim; ++j )
    {
        const Real delta = d.Get(j,0);
        if( delta >= Real(0) )
            d.Set(j,0,Real(1));
        else
            d.Set(j,0,Real(-1));
    }
    DiagonalScaleTrapezoid( RIGHT, UPPER, NORMAL, d, R, -iOff );
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
PanelHouseholder
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d )
{
    DEBUG_ONLY(
        CallStackEntry cse("rq::PanelHouseholder");
        if( A.Grid() != t.Grid() || t.Grid() != d.Grid() )
            LogicError("{A,t,d} must be distributed over the same grid");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int iOff = m-minDim;
    const Int jOff = n-minDim;

    t.SetRoot( A.DiagonalRoot(n-m) );
    d.SetRoot( A.DiagonalRoot(n-m) );
    t.AlignCols( A.DiagonalAlign(n-m) );
    d.AlignCols( A.DiagonalAlign(n-m) );
    t.Resize( minDim, 1 );
    d.Resize( minDim, 1 );

    const Grid& g = A.Grid();
    DistMatrix<F,STAR,MR  > a1L_STAR_MR(g);
    DistMatrix<F,MC,  STAR> z01_MC_STAR(g);

    for( Int k=minDim-1; k>=0; --k )
    {
        const Int ki = k + iOff;
        const Int kj = k + jOff;
        auto a10     = ViewRange( A, ki, 0,  ki+1, kj   );
        auto alpha11 = ViewRange( A, ki, kj, ki+1, kj+1 );
        auto A0L     = ViewRange( A, 0,  0,  ki,   kj+1 );
        auto a1L     = ViewRange( A, ki, 0,  ki+1, kj+1 );

        // Find tau and v such that
        //  |a10 alpha11| /I - tau |v^T| |conj(v) 1|\ = |0 beta|
        //                \        |1  |            /
        const F tau = RightReflector( alpha11, a10 );
        t.Set( k, 0, tau );

        // Temporarily set a1L = | v 1 |
        F alpha = 0;
        if( alpha11.IsLocal(0,0) )
        {
            alpha = alpha11.GetLocal(0,0);
            alpha11.SetLocal(0,0,1);
        }

        // A2R := A2R Hous(a1L^T,tau)
        //      = A2R (I - tau a1L^T conj(a1L))
        //      = A2R - tau (A2R a1L^T) conj(a1L)
        a1L_STAR_MR.AlignWith( A0L );
        a1L_STAR_MR = a1L;
        z01_MC_STAR.AlignWith( A0L );
        Zeros( z01_MC_STAR, A0L.Height(), 1 );
        LocalGemv( NORMAL, F(1), A0L, a1L_STAR_MR, F(0), z01_MC_STAR );
        z01_MC_STAR.SumOver( A0L.RowComm() );
        Ger
        ( -tau, z01_MC_STAR.LockedMatrix(), a1L_STAR_MR.LockedMatrix(),
          A0L.Matrix() ); 

        // Reset alpha11's value
        if( alpha11.IsLocal(0,0) )
            alpha11.SetLocal(0,0,alpha);
    }
    // Form d and rescale R
    auto R = View( A, 0, jOff, m, minDim );
    d = R.GetRealPartOfDiagonal();
    const Int diagLengthLoc = d.LocalHeight();
    typedef Base<F> Real;
    for( Int jLoc=0; jLoc<diagLengthLoc; ++jLoc )
    {
        const Real delta = d.GetLocal(jLoc,0);
        if( delta >= Real(0) )
            d.SetLocal(jLoc,0,Real(1));
        else
            d.SetLocal(jLoc,0,Real(-1));
    }
    DiagonalScaleTrapezoid( RIGHT, UPPER, NORMAL, d, R, -iOff );
}

template<typename F> 
inline void
PanelHouseholder( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("rq::PanelHouseholder"))
    DistMatrix<F,MD,STAR> t(A.Grid());
    DistMatrix<Base<F>,MD,STAR> d(A.Grid());
    PanelHouseholder( A, t, d );
}

} // namespace rq
} // namespace El

#endif // ifndef EL_RQ_PANEL_HPP
