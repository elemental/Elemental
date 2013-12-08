/*
   Copyright (c) 1992-2008 The University of Tennessee. 
   All rights reserved.

   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is partially based upon the LAPACK routines dlarfg.f and zlarfg.f.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_REFLECTOR_ROW_HPP
#define ELEM_LAPACK_REFLECTOR_ROW_HPP

#include "elemental/blas-like/level1/Nrm2.hpp"
#include "elemental/blas-like/level1/Scale.hpp"

namespace elem {
namespace reflector {

//
// Follows the LAPACK convention of defining tau such that
//
//   H = I - tau [1; w] [1, w'],
//
// but adjoint(H) [chi; y] = [beta; 0]. 
//
// In our case, y=x^T, and w=v^T, so that 
//
//   (I - conj(tau) [1; v^T] [1, conj(v)]) [chi; x^T] = [beta; 0],
//
// and thus
//
//   [chi, x] (I - conj(tau) [1; v^H] [1, v]) = [beta, 0].
//
// In the case of real data, everything simplifies to
//
//   [chi, x] (I - tau [1; v^T] [1, v]) = [beta, 0].
//
// On exit, chi is overwritten with beta, and x is overwritten with v.
//
// The major difference from LAPACK is in the treatment of the special case 
// of x=0, where LAPACK would put H := I, which is not a valid Householder 
// reflector. We instead follow the FLAME convention of defining H such that 
//    adjoint(H) [chi; 0] = [-chi; 0],
// which is accomplished by setting tau=2, and v^T=0.
//

template<typename R> 
inline R
Row( DistMatrix<R>& chi, DistMatrix<R>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("reflector::Row");
        if( chi.Grid() != x.Grid() )
            LogicError("chi and x must be distributed over the same grid");
        if( chi.Height() != 1 || chi.Width() != 1 )
            LogicError("chi must be a scalar");
        if( x.Height() != 1 )
            LogicError("x must be a row vector");
        if( chi.Grid().Row() != chi.ColAlign() )
            LogicError("Reflecting with incorrect row of processes");
        if( x.Grid().Row() != x.ColAlign() )
            LogicError("Reflecting with incorrect row of processes");
    )
    const Grid& grid = x.Grid();
    mpi::Comm rowComm = grid.RowComm();
    const Int gridCol = grid.Col();
    const Int gridWidth = grid.Width();
    const Int rowAlign = chi.RowAlign();

    std::vector<R> localNorms(gridWidth);
    R localNorm = Nrm2( x.LockedMatrix() ); 
    mpi::AllGather( &localNorm, 1, localNorms.data(), 1, rowComm );
    R norm = blas::Nrm2( gridWidth, localNorms.data(), 1 );

    if( norm == 0 )
    {
        if( gridCol == rowAlign )
            chi.SetLocal(0,0,-chi.GetLocal(0,0));
        return R(2);
    }

    R alpha;
    if( gridCol == rowAlign )
        alpha = chi.GetLocal(0,0);
    mpi::Broadcast( alpha, rowAlign, rowComm );

    R beta;
    if( alpha <= 0 )
        beta = lapack::SafeNorm( alpha, norm );
    else
        beta = -lapack::SafeNorm( alpha, norm );

    const R one = 1;
    const R safeMin = lapack::MachineSafeMin<R>();
    const R epsilon = lapack::MachineEpsilon<R>();
    const R safeInv = safeMin/epsilon;
    Int count = 0;
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

        localNorm = Nrm2( x.LockedMatrix() );
        mpi::AllGather( &localNorm, 1, localNorms.data(), 1, rowComm );
        norm = blas::Nrm2( gridWidth, localNorms.data(), 1 );
        if( alpha <= 0 )
            beta = lapack::SafeNorm( alpha, norm );
        else
            beta = -lapack::SafeNorm( alpha, norm );
    }

    R tau = (beta-alpha)/beta;
    Scale( one/(alpha-beta), x );

    for( Int j=0; j<count; ++j )
        beta *= safeInv;
    if( gridCol == rowAlign )
        chi.SetLocal(0,0,beta);
        
    return tau;
}

template<typename R>
inline Complex<R>
Row( DistMatrix<Complex<R> >& chi, DistMatrix<Complex<R> >& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("reflector::Row");    
        if( chi.Grid() != x.Grid() )
            LogicError("chi and x must be distributed over the same grid");
        if( chi.Height() != 1 || chi.Width() != 1 )
            LogicError("chi must be a scalar");
        if( x.Height() != 1 )
            LogicError("x must be a row vector");
        if( chi.Grid().Row() != chi.ColAlign() )
            LogicError("Reflecting with incorrect row of processes");
        if( x.Grid().Row() != x.ColAlign() )
            LogicError("Reflecting with incorrect row of processes");
    )
    typedef Complex<R> C;
    const Grid& grid = x.Grid();
    mpi::Comm rowComm = grid.RowComm();
    const Int gridCol = grid.Col();
    const Int gridWidth = grid.Width();
    const Int rowAlign = chi.RowAlign();

    std::vector<R> localNorms(gridWidth);
    R localNorm = Nrm2( x.LockedMatrix() ); 
    mpi::AllGather( &localNorm, 1, localNorms.data(), 1, rowComm );
    R norm = blas::Nrm2( gridWidth, localNorms.data(), 1 );

    C alpha;
    if( gridCol == rowAlign )
        alpha = chi.GetLocal(0,0);
    mpi::Broadcast( alpha, rowAlign, rowComm );

    if( norm == R(0) && alpha.imag() == R(0) )
    {
        if( gridCol == rowAlign )
            chi.SetLocal(0,0,-chi.GetLocal(0,0));
        return C(2);
    }

    R beta;
    if( alpha.real() <= 0 )
        beta = lapack::SafeNorm( alpha.real(), alpha.imag(), norm );
    else
        beta = -lapack::SafeNorm( alpha.real(), alpha.imag(), norm );

    const R one = 1;
    const R safeMin = lapack::MachineSafeMin<R>();
    const R epsilon = lapack::MachineEpsilon<R>();
    const R safeInv = safeMin/epsilon;
    Int count = 0;
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

        localNorm = Nrm2( x.LockedMatrix() );
        mpi::AllGather( &localNorm, 1, localNorms.data(), 1, rowComm );
        norm = blas::Nrm2( gridWidth, localNorms.data(), 1 );
        if( alpha.real() <= 0 )
            beta = lapack::SafeNorm( alpha.real(), alpha.imag(), norm );
        else
            beta = -lapack::SafeNorm( alpha.real(), alpha.imag(), norm );
    }

    C tau = C( (beta-alpha.real())/beta, -alpha.imag()/beta );
    Scale( one/(alpha-beta), x );

    for( Int j=0; j<count; ++j )
        beta *= safeInv;
    if( gridCol == rowAlign )
        chi.SetLocal(0,0,beta);
        
    return tau;
}

} // namespace reflector
} // namespace elem

#endif // ifndef ELEM_LAPACK_REFLECTOR_ROW_HPP
