/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LQ_PANELHOUSEHOLDER_HPP
#define EL_LQ_PANELHOUSEHOLDER_HPP



namespace El {
namespace lq {

template<typename F>
inline void
PanelHouseholder( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d )
{
    DEBUG_ONLY(CallStackEntry cse("lq::PanelHouseholder"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    t.Resize( minDim, 1 );
    d.Resize( minDim, 1 );

    Matrix<F> z21;

    for( Int k=0; k<minDim; ++k )
    {
        auto alpha11 = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12     = ViewRange( A, k,   k+1, k+1, n   );
        auto a1R     = ViewRange( A, k,   k,   k+1, n   );
        auto A2R     = ViewRange( A, k+1, k,   m,   n   );

        // Find tau and v such that
        //  |alpha11 a12| /I - tau |1  | |1 conj(v)|\ = |beta 0|
        //                \        |v^T|            /
        const F tau = RightReflector( alpha11, a12 );
        t.Set( k, 0, tau );

        // Temporarily set a1R = | 1 v |
        const F alpha = alpha11.Get(0,0);
        alpha11.Set(0,0,1);

        // A2R := A2R Hous(a1R^T,tau)
        //      = A2R (I - tau a1R^T conj(a1R))
        //      = A2R - tau (A2R a1R^T) conj(a1R)
        Zeros( z21, A2R.Height(), 1 );
        Gemv( NORMAL, F(1), A2R, a1R, F(0), z21 );
        Ger( -tau, z21, a1R, A2R );

        // Reset alpha11's value
        alpha11.Set(0,0,alpha);
    }
    // Form d and rescale L
    auto L = View( A, 0, 0, m, minDim );
    d = L.GetRealPartOfDiagonal();
    typedef Base<F> Real;
    for( Int j=0; j<minDim; ++j )
    {
        const Real delta = d.Get(j,0);
        if( delta >= Real(0) )
            d.Set(j,0,Real(1));
        else
            d.Set(j,0,Real(-1));
    }
    DiagonalScaleTrapezoid( RIGHT, LOWER, NORMAL, d, L );
}

template<typename F>
inline void
PanelHouseholder( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("lq::PanelHouseholder"))
    Matrix<F> t;
    Matrix<Base<F>> d;
    PanelHouseholder( A, t, d );
}

template<typename F>
inline void
PanelHouseholder
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d )
{
    DEBUG_ONLY(
        CallStackEntry cse("lq::PanelHouseholder");
        if( A.Grid() != t.Grid() || t.Grid() != d.Grid() )
            LogicError("{A,t,d} must be distributed over the same grid");
    )
    t.SetRoot( A.DiagonalRoot() );
    d.SetRoot( A.DiagonalRoot() );
    t.AlignCols( A.DiagonalAlign() );
    d.AlignCols( A.DiagonalAlign() );
    const Grid& g = A.Grid();
    DistMatrix<F,STAR,MR  > a1R_STAR_MR(g);
    DistMatrix<F,MC,  STAR> z21_MC_STAR(g);

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    t.Resize( minDim, 1 );
    d.Resize( minDim, 1 );

    for( Int k=0; k<minDim; ++k )
    {
        auto alpha11 = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12     = ViewRange( A, k,   k+1, k+1, n   );
        auto a1R     = ViewRange( A, k,   k,   k+1, n   );
        auto A2R     = ViewRange( A, k+1, k,   m,   n   );

        // Find tau and v such that
        //  |alpha11 a12| /I - tau |1  | |1 conj(v)|\ = |beta 0|
        //                \        |v^T|            /
        const F tau = RightReflector( alpha11, a12 );
        t.Set( k, 0, tau );

        // Temporarily set a1R = | 1 v |
        F alpha = 0;
        if( alpha11.IsLocal(0,0) )
        {
            alpha = alpha11.GetLocal(0,0);
            alpha11.SetLocal(0,0,1);
        }

        // A2R := A2R Hous(a1R^T,tau)
        //      = A2R (I - tau a1R^T conj(a1R))
        //      = A2R - tau (A2R a1R^T) conj(a1R)
        a1R_STAR_MR.AlignWith( A2R );
        a1R_STAR_MR = a1R;
        z21_MC_STAR.AlignWith( A2R );
        Zeros( z21_MC_STAR, A2R.Height(), 1 );
        LocalGemv( NORMAL, F(1), A2R, a1R_STAR_MR, F(0), z21_MC_STAR );
        z21_MC_STAR.SumOver( A2R.RowComm() );
        Ger
        ( -tau, z21_MC_STAR.LockedMatrix(), a1R_STAR_MR.LockedMatrix(),
          A2R.Matrix() );

        // Reset alpha11's value
        if( alpha11.IsLocal(0,0) )
            alpha11.SetLocal(0,0,alpha);
    }
    // Form d and rescale L
    auto L = View( A, 0, 0, m, minDim );
    d = L.GetRealPartOfDiagonal();
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
    DiagonalScaleTrapezoid( RIGHT, LOWER, NORMAL, d, L );
}

template<typename F>
inline void
PanelHouseholder( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("lq::PanelHouseholder"))
    DistMatrix<F,MD,STAR> t(A.Grid());
    DistMatrix<Base<F>,MD,STAR> d(A.Grid());
    PanelHouseholder( A, t, d );
}

} // namespace lq
} // namespace El

#endif // ifndef EL_LQ_PANELHOUSEHOLDER_HPP
