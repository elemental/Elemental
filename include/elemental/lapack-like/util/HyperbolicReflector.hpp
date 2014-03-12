/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HYPERBOLICREFLECTOR_HPP
#define ELEM_HYPERBOLICREFLECTOR_HPP

#include ELEM_CONJUGATE_INC
#include ELEM_NRM2_INC
#include ELEM_SCALE_INC
#include ELEM_ZERO_INC

#include "./HyperbolicReflector/Col.hpp"
#include "./HyperbolicReflector/Row.hpp"

namespace elem {

//
// (I - 1/tau w w^H Sigma) x = -lambda e_0
// 
// where Sigma is diag(+1,-1,...,-1)

template<typename F>
inline F
LeftHyperbolicReflector( Matrix<F>& chi, Matrix<F>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("LeftHyperbolicReflector");
        if( chi.Height() != 1 || chi.Width() != 1 )
            LogicError("chi must be a scalar");
        if( x.Height() != 1 && x.Width() != 1 )
            LogicError("x must be a vector");
        if( ImagPart(chi.Get(0,0)) != Base<F>(0) )
            LogicError("chi is assumed to be real");
    )

    // Compute lambda = sgn(chi) sqrt([chi;x]^H Sigma [chi;x])
    //                = sgn(chi) sqrt(chi^2 - x^H x)    
    typedef Base<F> Real;
    const Real alpha = chi.GetRealPart(0,0);
    const Real xNrm = Nrm2( x ); 
    const Real delta = alpha*alpha - xNrm*xNrm;
    if( delta < Real(0) )
        LogicError("Attempted to square-root a negative number");
    const Real lambda = ( alpha>=0 ? Sqrt(delta) : -Sqrt(delta) );
    chi.Set(0,0,-lambda);

    // Implicitly define 
    //     w := [chi;x] + lambda e_0, and 
    //     kappa = chi + lambda,
    // so that
    //     tau := (w^H Sigma w) / 2 = (delta + lambda^2 + 2 chi lambda) / 2
    //          = delta + chi lambda
    // then normalize w so that its first entry is one
    // TODO: Introduce a threshold instead of the approach from 
    //       van de Geijn and van Zee's "High-performance up-and-downdating
    //       via Householder-like transformations"
    const Real kappa = alpha + lambda;
    if( kappa == Real(0) )
    {
        Zero( x );        
        return Real(1);
    }
    else
    {
        Scale( Real(1)/kappa, x );
        return (delta+alpha*lambda)/(kappa*kappa);
    }
}

template<typename F>
inline F
LeftHyperbolicReflector( F& chi, Matrix<F>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("LeftHyperbolicReflector");
        if( x.Height() != 1 && x.Width() != 1 )
            LogicError("x must be a vector");
        if( ImagPart(chi) != Base<F>(0) )
            LogicError("chi is assumed to be real");
    )

    // Compute lambda = sgn(chi) sqrt([chi;x]^H Sigma [chi;x])
    //                = sgn(chi) sqrt(chi^2 - x^H x)    
    typedef Base<F> Real;
    const Real alpha = RealPart(chi);
    const Real xNrm = Nrm2( x ); 
    const Real delta = alpha*alpha - xNrm*xNrm;
    if( delta < Real(0) )
        LogicError("Attempted to square-root a negative number");
    const Real lambda = ( alpha>=0 ? Sqrt(delta) : -Sqrt(delta) );
    chi = -lambda;

    // Implicitly define 
    //     w := [chi;x] + lambda e_0, and 
    //     kappa = chi + lambda,
    // so that
    //     tau := (w^H Sigma w) / 2 = (delta + lambda^2 + 2 chi lambda) / 2
    //          = delta + chi lambda
    // then normalize w so that its first entry is one
    // TODO: Introduce a threshold instead of the approach from 
    //       van de Geijn and van Zee's "High-performance up-and-downdating
    //       via Householder-like transformations"
    const Real kappa = alpha + lambda;
    if( kappa == Real(0) )
    {
        Zero( x );        
        return Real(1);
    }
    else
    {
        Scale( Real(1)/kappa, x );
        return (delta+alpha*lambda)/(kappa*kappa);
    }
}

template<typename F,Dist U,Dist V>
inline F
LeftHyperbolicReflector( DistMatrix<F,U,V>& chi, DistMatrix<F,U,V>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("LeftHyperbolicReflector");
        if( chi.Grid() != x.Grid() )
            LogicError("chi and x must be distributed over the same grid");
        if( chi.Height() != 1 || chi.Width() != 1 )
            LogicError("chi must be a scalar");
        if( x.Width() != 1 )
            LogicError("x must be a column vector");
    )
    F tau;
    if( x.RowRank() == x.RowAlign() )
        tau = hyp_reflector::Col( chi, x );
    mpi::Broadcast( tau, x.RowAlign(), x.RowComm() );
    return tau;
}

template<typename F,Dist U,Dist V>
inline F
LeftHyperbolicReflector( F& chi, DistMatrix<F,U,V>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("LeftHyperbolicReflector");
        if( x.Width() != 1 )
            LogicError("x must be a column vector");
    )
    F tau;
    if( x.RowRank() == x.RowAlign() )
        tau = hyp_reflector::Col( chi, x );
    mpi::Broadcast( tau, x.RowAlign(), x.RowComm() );
    return tau;
}

//
// x (I - 1/tau Sigma w^H w) = -lambda e_0^T
// 
// where Sigma is diag(+1,-1,...,-1)

template<typename F>
inline F
RightHyperbolicReflector( Matrix<F>& chi, Matrix<F>& x )
{
    DEBUG_ONLY(CallStackEntry cse("RightHyperbolicReflector"))
    const F tau = LeftHyperbolicReflector( chi, x );
    Conjugate( x );
    return tau;
}

template<typename F>
inline F
RightHyperbolicReflector( F& chi, Matrix<F>& x )
{
    DEBUG_ONLY(CallStackEntry cse("RightHyperbolicReflector"))
    const F tau = LeftHyperbolicReflector( chi, x );
    Conjugate( x );
    return tau;
}

template<typename F,Dist U,Dist V>
inline F
RightHyperbolicReflector( DistMatrix<F,U,V>& chi, DistMatrix<F,U,V>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("RightHyperbolicReflector");
        if( chi.Grid() != x.Grid() )
            LogicError("chi and x must be distributed over the same grid");
        if( chi.Height() != 1 || chi.Width() != 1 )
            LogicError("chi must be a scalar");
        if( x.Height() != 1 )
            LogicError("x must be a row vector");
    )
    F tau;
    if( x.ColRank() == x.ColAlign() )
        tau = hyp_reflector::Row( chi, x );
    mpi::Broadcast( tau, x.ColAlign(), x.ColComm() );
    return tau;
}

template<typename F,Dist U,Dist V>
inline F
RightHyperbolicReflector( F& chi, DistMatrix<F,U,V>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("RightHyperbolicReflector");
        if( x.Height() != 1 )
            LogicError("x must be a row vector");
    )
    F tau;
    if( x.ColRank() == x.ColAlign() )
        tau = hyp_reflector::Row( chi, x );
    mpi::Broadcast( tau, x.ColAlign(), x.ColComm() );
    return tau;
}

} // namespace elem

#endif // ifndef ELEM_HYPERBOLICREFLECTOR_HPP
