/*
   Copyright (C) 1992-2008 The University of Tennessee
   All rights reserved.

   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is partially based upon the LAPACK routines dlarfg.f and zlarfg.f.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_REFLECTOR_HPP
#define LAPACK_REFLECTOR_HPP

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

template<typename R>
inline R
Reflector( Matrix<R>& chi, Matrix<R>& x )
{
#ifndef RELEASE
    PushCallStack("Reflector");
#endif
    R norm = Nrm2( x );
    if( norm == 0 )
    {
        chi.Set(0,0,-chi.Get(0,0));
#ifndef RELEASE
        PopCallStack();
#endif
        return R(2);
    }

    R beta;
    R alpha = chi.Get(0,0);
    if( alpha <= 0 )
        beta = lapack::SafeNorm( alpha, norm );
    else
        beta = -lapack::SafeNorm( alpha, norm );

    const R one = 1;
    const R safeMin = lapack::MachineSafeMin<R>();
    const R epsilon = lapack::MachineEpsilon<R>();
    const R safeInv = safeMin/epsilon;
    int count = 0;
    if( Abs(beta) < safeInv )
    {
        R invOfSafeInv = one/safeInv;
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

    R tau = (beta-alpha) / beta;
    Scale( one/(alpha-beta), x );

    for( int j=0; j<count; ++j )
        beta *= safeInv;
    chi.Set(0,0,beta);
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}

template<typename R>
inline R
Reflector( R& chi, int m, R* x, int incx )
{
#ifndef RELEASE
    PushCallStack("Reflector");
#endif
    R norm = blas::Nrm2( m, x, incx );
    if( norm == 0 )
    {
        x[0] = -x[0];
#ifndef RELEASE
        PopCallStack();
#endif
        return R(2);
    }

    R beta;
    R alpha = chi;
    if( alpha <= 0 )
        beta = lapack::SafeNorm( alpha, norm );
    else
        beta = -lapack::SafeNorm( alpha, norm );

    // Avoid overflow by scaling the vector
    const R one = 1;
    const R safeMin = lapack::MachineSafeMin<R>();
    const R epsilon = lapack::MachineEpsilon<R>();
    const R safeInv = safeMin/epsilon;
    int count=0;
    if( Abs(beta) < safeInv )
    {
        R invOfSafeInv = one/safeInv;
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

    R tau = (beta-alpha) / beta;
    blas::Scal( m, one/(alpha-beta), x, incx );

    for( int j=0; j<count; ++j )
        beta *= safeInv;
    chi = beta;
#ifndef RELEASE
    PopCallStack();
#endif
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

template<typename R>
inline Complex<R>
Reflector( Matrix<Complex<R> >& chi, Matrix<Complex<R> >& x )
{
#ifndef RELEASE
    PushCallStack("Reflector");
#endif
    typedef Complex<R> C;

    R norm = Nrm2( x );
    C alpha = chi.Get(0,0);

    if( norm == 0 && alpha.imag == R(0) )
    {
        chi.Set(0,0,-chi.Get(0,0));
#ifndef RELEASE
        PopCallStack();
#endif
        return C(2);
    }

    R beta;
    if( alpha.real <= 0 )
        beta = lapack::SafeNorm( alpha.real, alpha.imag, norm );
    else
        beta = -lapack::SafeNorm( alpha.real, alpha.imag, norm );

    const R one = 1;
    const R safeMin = lapack::MachineSafeMin<R>();
    const R epsilon = lapack::MachineEpsilon<R>();
    const R safeInv = safeMin/epsilon;
    int count = 0;
    if( Abs(beta) < safeInv )
    {
        R invOfSafeInv = one/safeInv;
        do
        {
            ++count;
            Scale( invOfSafeInv, x );
            alpha *= invOfSafeInv;
            beta *= invOfSafeInv;
        } while( Abs(beta) < safeInv );

        norm = Nrm2( x );
        if( alpha.real <= 0 )
            beta = lapack::SafeNorm( alpha.real, alpha.imag, norm );
        else
            beta = -lapack::SafeNorm( alpha.real, alpha.imag, norm );
    }

    C tau = C( (beta-alpha.real)/beta, -alpha.imag/beta );
    Scale( one/(alpha-beta), x );

    for( int j=0; j<count; ++j )
        beta *= safeInv;
    chi.Set(0,0,beta);
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}

template<typename R>
inline Complex<R>
Reflector( Complex<R>& chi, int m, Complex<R>* x, int incx )
{
#ifndef RELEASE
    PushCallStack("Reflector");
#endif
    typedef Complex<R> C;

    R norm = blas::Nrm2( m, x, incx );
    C alpha = chi;

    if( norm == 0 && alpha.imag == R(0) )
    {
        chi = -chi;
#ifndef RELEASE
        PopCallStack();
#endif
        return C(2);
    }

    R beta;
    if( alpha.real <= 0 )
        beta = lapack::SafeNorm( alpha.real, alpha.imag, norm );
    else
        beta = -lapack::SafeNorm( alpha.real, alpha.imag, norm );

    const R one = 1;
    const R safeMin = lapack::MachineSafeMin<R>();
    const R epsilon = lapack::MachineEpsilon<R>();
    const R safeInv = safeMin/epsilon;
    int count = 0;
    if( Abs(beta) < safeInv )
    {
        R invOfSafeInv = one/safeInv;
        do
        {
            ++count;
            blas::Scal( m, invOfSafeInv, x, incx );
            alpha *= invOfSafeInv;
            beta *= invOfSafeInv;
        } while( Abs(beta) < safeInv );

        norm = blas::Nrm2( m, x, incx );
        if( alpha.real <= 0 )
            beta = lapack::SafeNorm( alpha.real, alpha.imag, norm );
        else
            beta = -lapack::SafeNorm( alpha.real, alpha.imag, norm );
    }

    C tau = C( (beta-alpha.real)/beta, -alpha.imag/beta );
    blas::Scal( m, one/(alpha-beta), x, incx );

    for( int j=0; j<count; ++j )
        beta *= safeInv;
    chi = beta;

#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}

template<typename F>
inline F
Reflector( DistMatrix<F>& chi, DistMatrix<F>& x )
{
#ifndef RELEASE
    PushCallStack("Reflector");
    if( chi.Grid() != x.Grid() )
        throw std::logic_error
        ("chi and x must be distributed over the same grid");
    if( chi.Height() != 1 || chi.Width() != 1 )
        throw std::logic_error("chi must be a scalar");
    if( x.Height() != 1 && x.Width() != 1 )
        throw std::logic_error("x must be a vector");
#endif
    const Grid& g = x.Grid();
    F tau;
    if( x.Width() == 1 && x.RowAlignment() == chi.RowAlignment() )
    {
        if( g.Col() == x.RowAlignment() )
            tau = reflector::Col( chi, x );
        mpi::Broadcast( &tau, 1, x.RowAlignment(), g.RowComm() );
    }
    else
    {
        if( g.Row() == x.ColAlignment() )
            tau = reflector::Row( chi, x );
        mpi::Broadcast( &tau, 1, x.ColAlignment(), g.ColComm() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}

} // namespace elem

#endif // ifndef LAPACK_REFLECTOR_HPP
