/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_QR_PANEL_HPP
#define EL_QR_PANEL_HPP

#include EL_ZEROS_INC

namespace El {
namespace qr {

template<typename F> 
inline void
PanelHouseholder( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d )
{
    DEBUG_ONLY(CallStackEntry cse("qr::PanelHouseholder"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    t.Resize( minDim, 1 );
    d.Resize( minDim, 1 );

    Matrix<F> z21;

    for( Int k=0; k<minDim; ++k )
    {
        auto alpha11 = ViewRange( A, k,   k,   k+1, k+1 );
        auto a21     = ViewRange( A, k+1, k,   m,   k+1 );
        auto aB1     = ViewRange( A, k,   k,   m,   k+1 );
        auto AB2     = ViewRange( A, k,   k+1, m,   n   );

        // Find tau and u such that
        //  / I - tau | 1 | | 1, u^H | \ | alpha11 | = | beta |
        //  \         | u |            / |     a21 | = |    0 |
        const F tau = LeftReflector( alpha11, a21 );
        t.Set( k, 0, tau );

        // Temporarily set aB1 = | 1 |
        //                       | u |
        const F alpha = alpha11.Get(0,0);
        alpha11.Set(0,0,1);

        // AB2 := Hous(aB1,tau) AB2
        //      = (I - tau aB1 aB1^H) AB2
        //      = AB2 - tau aB1 (AB2^H aB1)^H
        Zeros( z21, AB2.Width(), 1 );
        Gemv( ADJOINT, F(1), AB2, aB1, F(0), z21 );
        Ger( -tau, aB1, z21, AB2 );

        // Replace alpha11's value
        alpha11.Set(0,0,alpha);
    }
    // Form d and rescale R
    auto R = View( A, 0, 0, minDim, n );
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
    DiagonalScaleTrapezoid( LEFT, UPPER, NORMAL, d, R );
}

template<typename F> 
inline void
PanelHouseholder( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("qr::PanelHouseholder"))
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
        CallStackEntry cse("qr::PanelHouseholder");
        if( A.Grid() != t.Grid() || t.Grid() != d.Grid() )
            LogicError("{A,t,d} must be distributed over the same grid");
    )
    t.SetRoot( A.DiagonalRoot() );
    d.SetRoot( A.DiagonalRoot() );
    t.AlignCols( A.DiagonalAlign() );
    d.AlignCols( A.DiagonalAlign() );
    const Grid& g = A.Grid();
    DistMatrix<F,MC,STAR> aB1_MC_STAR(g);
    DistMatrix<F,MR,STAR> z21_MR_STAR(g);

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    t.Resize( minDim, 1 );
    d.Resize( minDim, 1 );

    for( Int k=0; k<minDim; ++k )
    {
        auto alpha11 = ViewRange( A, k,   k,   k+1, k+1 );
        auto a21     = ViewRange( A, k+1, k,   m,   k+1 );
        auto aB1     = ViewRange( A, k,   k,   m,   k+1 );
        auto AB2     = ViewRange( A, k,   k+1, m,   n   );

        // Find tau and u such that
        //  / I - tau | 1 | | 1, u^H | \ | alpha11 | = | beta |
        //  \         | u |            / |     a21 | = |    0 |
        const F tau = LeftReflector( alpha11, a21 );
        t.Set( k, 0, tau );

        // Temporarily set aB1 = | 1 |
        //                       | u |
        F alpha = 0;
        if( alpha11.IsLocal(0,0) )
        {
            alpha = alpha11.GetLocal(0,0);
            alpha11.SetLocal(0,0,F(1));
        }

        // AB2 := Hous(aB1,tau) AB2
        //      = (I - tau aB1 aB1^H) AB2
        //      = AB2 - tau aB1 (AB2^H aB1)^H
        aB1_MC_STAR.AlignWith( AB2 );
        aB1_MC_STAR = aB1;
        z21_MR_STAR.AlignWith( AB2 );
        Zeros( z21_MR_STAR, AB2.Width(), 1 );
        LocalGemv( ADJOINT, F(1), AB2, aB1_MC_STAR, F(0), z21_MR_STAR );
        z21_MR_STAR.SumOver( AB2.ColComm() );
        Ger
        ( -tau, aB1_MC_STAR.LockedMatrix(), z21_MR_STAR.LockedMatrix(),
          AB2.Matrix() );

        // Replace alpha11's value
        if( alpha11.IsLocal(0,0) )
            alpha11.SetLocal(0,0,alpha);
    }
    // Form d and rescale R
    auto R = View( A, 0, 0, minDim, n );
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
    DiagonalScaleTrapezoid( LEFT, UPPER, NORMAL, d, R );
}

template<typename F> 
inline void
PanelHouseholder( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("qr::PanelHouseholder"))
    DistMatrix<F,MD,STAR> t(A.Grid());
    PanelHouseholder( A, t );
}

} // namespace qr
} // namespace El

#endif // ifndef EL_QR_PANEL_HPP
