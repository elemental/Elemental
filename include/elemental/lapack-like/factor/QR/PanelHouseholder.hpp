/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_QR_PANEL_HPP
#define ELEM_QR_PANEL_HPP

#include ELEM_GEMV_INC
#include ELEM_GER_INC
#include ELEM_REFLECTOR_INC
#include ELEM_ZEROS_INC

namespace elem {
namespace qr {

template<typename F> 
inline void
PanelHouseholder( Matrix<F>& A, Matrix<F>& t )
{
    DEBUG_ONLY(CallStackEntry cse("qr::PanelHouseholder"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    t.Resize( minDim, 1 );

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
}

template<typename F> 
inline void
PanelHouseholder( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("qr::PanelHouseholder"))
    Matrix<F> t;
    PanelHouseholder( A, t );
}

template<typename F> 
inline void
PanelHouseholder( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t )
{
    DEBUG_ONLY(
        CallStackEntry cse("qr::PanelHouseholder");
        if( A.Grid() != t.Grid() )
            LogicError("{A,t} must be distributed over the same grid");
        if( !A.DiagonalAlignedWith( t, 0 ) )
            LogicError("t must be aligned with A's main diagonal");
    )
    const Grid& g = A.Grid();
    DistMatrix<F,MC,STAR> aB1_MC_STAR(g);
    DistMatrix<F,MR,STAR> z21_MR_STAR(g);

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    t.Resize( minDim, 1 );

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
} // namespace elem

#endif // ifndef ELEM_QR_PANEL_HPP
