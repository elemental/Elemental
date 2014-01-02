/*
   Copyright (C) 1992-2008 The University of Tennessee
   All rights reserved.

   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is partially based upon the LAPACK routines dlarfg.f and zlarfg.f.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_REFLECTOR_HPP
#define ELEM_LAPACK_REFLECTOR_HPP

#include "elemental/blas-like/level1/Nrm2.hpp"
#include "elemental/blas-like/level1/Scale.hpp"

#include "./Reflector/Col.hpp"
#include "./Reflector/Row.hpp"

namespace elem {

//
// Follows the LAPACK convention of defining tau such that
//
//   H = I - tau [1; v] [1, v'],
//
// but adjoint(H) [chi; x] = [beta; 0]. 
//
// Note that the adjoint of H is applied. In this case, where the data is real,
// H' = H, so there is no complication.
//
// On exit, chi is overwritten with beta, and x is overwritten with v.
//
// The major difference from LAPACK is in the treatment of the special case 
// of x=0, where LAPACK would put H := I, which is not a valid Householder 
// reflector. We instead follow the FLAME convention of defining H such that 
//    adjoint(H) [chi; 0] = [-chi; 0],
// which is accomplished by setting tau=2, and v=0.
//

template<typename Real>
inline Real
Reflector( Matrix<Real>& chi, Matrix<Real>& x )
{
    DEBUG_ONLY(CallStackEntry cse("Reflector"))
    Real norm = Nrm2( x );
    if( norm == Real(0) )
    {
        chi.Set(0,0,-chi.Get(0,0));
        return Real(2);
    }

    Real beta;
    Real alpha = chi.Get(0,0);
    if( alpha <= 0 )
        beta = lapack::SafeNorm( alpha, norm );
    else
        beta = -lapack::SafeNorm( alpha, norm );

    const Real one = 1;
    const Real safeMin = lapack::MachineSafeMin<Real>();
    const Real epsilon = lapack::MachineEpsilon<Real>();
    const Real safeInv = safeMin/epsilon;
    Int count = 0;
    if( Abs(beta) < safeInv )
    {
        Real invOfSafeInv = one/safeInv;
        do
        {
            ++count;
            Scale( invOfSafeInv, x );
            alpha *= invOfSafeInv;
            beta *= invOfSafeInv;
        } while( Abs(beta) < safeInv );

        norm = Nrm2( x );
        if( alpha <= 0 )
            beta = lapack::SafeNorm( alpha, norm );
        else
            beta = -lapack::SafeNorm( alpha, norm );
    }

    Real tau = (beta-alpha) / beta;
    Scale( one/(alpha-beta), x );

    for( Int j=0; j<count; ++j )
        beta *= safeInv;
    chi.Set(0,0,beta);
    return tau;
}

template<typename Real>
inline Real
Reflector( Real& chi, Int m, Real* x, Int incx )
{
    DEBUG_ONLY(CallStackEntry cse("Reflector"))
    Real norm = blas::Nrm2( m, x, incx );
    if( norm == 0 )
    {
        x[0] = -x[0];
        return Real(2);
    }

    Real beta;
    Real alpha = chi;
    if( alpha <= 0 )
        beta = lapack::SafeNorm( alpha, norm );
    else
        beta = -lapack::SafeNorm( alpha, norm );

    // Avoid overflow by scaling the vector
    const Real one = 1;
    const Real safeMin = lapack::MachineSafeMin<Real>();
    const Real epsilon = lapack::MachineEpsilon<Real>();
    const Real safeInv = safeMin/epsilon;
    Int count=0;
    if( Abs(beta) < safeInv )
    {
        Real invOfSafeInv = one/safeInv;
        do
        {
            ++count;
            blas::Scal( m, invOfSafeInv, x, incx );
            alpha *= invOfSafeInv;
            beta *= invOfSafeInv;
        } while( Abs(beta) < safeInv );

        norm = blas::Nrm2( m, x, incx );
        if( alpha <= 0 )
            beta = lapack::SafeNorm( alpha, norm );
        else
            beta = -lapack::SafeNorm( alpha, norm );
    }

    Real tau = (beta-alpha) / beta;
    blas::Scal( m, one/(alpha-beta), x, incx );

    for( Int j=0; j<count; ++j )
        beta *= safeInv;
    chi = beta;
    return tau;
}

//
// Follows the LAPACK convention of defining tau such that
//
//   H = I - tau [1; v] [1, v'],
//
// but adjoint(H) [chi; x] = [beta; 0]. 
//
// Note that the adjoint of H is applied. 
//
// On exit, chi is overwritten with beta, and x is overwritten with v.
//
// The major difference from LAPACK is in the treatment of the special case 
// of x=0, where LAPACK would put H := I, which is not a valid Householder 
// reflector. We instead follow the FLAME convention of defining H such that 
//    adjoint(H) [chi; 0] = [-chi; 0],
// which is accomplished by setting tau=2, and v=0.
//

template<typename Real>
inline Complex<Real>
Reflector( Matrix<Complex<Real> >& chi, Matrix<Complex<Real> >& x )
{
    DEBUG_ONLY(CallStackEntry cse("Reflector"))
    typedef Complex<Real> C;

    Real norm = Nrm2( x );
    C alpha = chi.Get(0,0);

    if( norm == Real(0) && alpha.imag() == Real(0) )
    {
        chi.Set(0,0,-chi.Get(0,0));
        return C(2);
    }

    Real beta;
    if( alpha.real() <= 0 )
        beta = lapack::SafeNorm( alpha.real(), alpha.imag(), norm );
    else
        beta = -lapack::SafeNorm( alpha.real(), alpha.imag(), norm );

    const Real one = 1;
    const Real safeMin = lapack::MachineSafeMin<Real>();
    const Real epsilon = lapack::MachineEpsilon<Real>();
    const Real safeInv = safeMin/epsilon;
    Int count = 0;
    if( Abs(beta) < safeInv )
    {
        Real invOfSafeInv = one/safeInv;
        do
        {
            ++count;
            Scale( invOfSafeInv, x );
            alpha *= invOfSafeInv;
            beta *= invOfSafeInv;
        } while( Abs(beta) < safeInv );

        norm = Nrm2( x );
        if( alpha.real() <= 0 )
            beta = lapack::SafeNorm( alpha.real(), alpha.imag(), norm );
        else
            beta = -lapack::SafeNorm( alpha.real(), alpha.imag(), norm );
    }

    C tau = C( (beta-alpha.real())/beta, -alpha.imag()/beta );
    Scale( one/(alpha-beta), x );

    for( Int j=0; j<count; ++j )
        beta *= safeInv;
    chi.Set(0,0,beta);
    return tau;
}

template<typename Real>
inline Complex<Real>
Reflector( Complex<Real>& chi, Int m, Complex<Real>* x, Int incx )
{
    DEBUG_ONLY(CallStackEntry cse("Reflector"))
    typedef Complex<Real> C;

    Real norm = blas::Nrm2( m, x, incx );
    C alpha = chi;

    if( norm == Real(0) && alpha.imag() == Real(0) )
    {
        chi = -chi;
        return C(2);
    }

    Real beta;
    if( alpha.real() <= 0 )
        beta = lapack::SafeNorm( alpha.real(), alpha.imag(), norm );
    else
        beta = -lapack::SafeNorm( alpha.real(), alpha.imag(), norm );

    const Real one = 1;
    const Real safeMin = lapack::MachineSafeMin<Real>();
    const Real epsilon = lapack::MachineEpsilon<Real>();
    const Real safeInv = safeMin/epsilon;
    Int count = 0;
    if( Abs(beta) < safeInv )
    {
        Real invOfSafeInv = one/safeInv;
        do
        {
            ++count;
            blas::Scal( m, invOfSafeInv, x, incx );
            alpha *= invOfSafeInv;
            beta *= invOfSafeInv;
        } while( Abs(beta) < safeInv );

        norm = blas::Nrm2( m, x, incx );
        if( alpha.real() <= 0 )
            beta = lapack::SafeNorm( alpha.real(), alpha.imag(), norm );
        else
            beta = -lapack::SafeNorm( alpha.real(), alpha.imag(), norm );
    }

    C tau = C( (beta-alpha.real())/beta, -alpha.imag()/beta );
    blas::Scal( m, one/(alpha-beta), x, incx );

    for( Int j=0; j<count; ++j )
        beta *= safeInv;
    chi = beta;

    return tau;
}

template<typename F>
inline F
Reflector( DistMatrix<F>& chi, DistMatrix<F>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("Reflector");
        if( chi.Grid() != x.Grid() )
            LogicError("chi and x must be distributed over the same grid");
        if( chi.Height() != 1 || chi.Width() != 1 )
            LogicError("chi must be a scalar");
        if( x.Height() != 1 && x.Width() != 1 )
            LogicError("x must be a vector");
    )
    const Grid& g = x.Grid();
    F tau;
    if( x.Width() == 1 && x.RowAlign() == chi.RowAlign() )
    {
        if( g.Col() == x.RowAlign() )
            tau = reflector::Col( chi, x );
        mpi::Broadcast( tau, x.RowAlign(), g.RowComm() );
    }
    else
    {
        if( g.Row() == x.ColAlign() )
            tau = reflector::Row( chi, x );
        mpi::Broadcast( tau, x.ColAlign(), g.ColComm() );
    }
    return tau;
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_REFLECTOR_HPP
