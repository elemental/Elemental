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
#ifndef LAPACK_REFLECTOR_ROW_HPP
#define LAPACK_REFLECTOR_ROW_HPP

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
#ifndef RELEASE
    PushCallStack("reflector::Row");
    if( chi.Grid() != x.Grid() )
        throw std::logic_error
        ("chi and x must be distributed over the same grid");
    if( chi.Height() != 1 || chi.Width() != 1 )
        throw std::logic_error("chi must be a scalar");
    if( x.Height() != 1 )
        throw std::logic_error("x must be a row vector");
    if( chi.Grid().Row() != chi.ColAlignment() )
        throw std::logic_error("Reflecting with incorrect row of processes");
    if( x.Grid().Row() != x.ColAlignment() )
        throw std::logic_error("Reflecting with incorrect row of processes");
#endif
    const Grid& grid = x.Grid();
    mpi::Comm rowComm = grid.RowComm();
    const int gridCol = grid.Col();
    const int gridWidth = grid.Width();
    const int rowAlignment = chi.RowAlignment();

    std::vector<R> localNorms(gridWidth);
    R localNorm = Nrm2( x.LockedMatrix() ); 
    mpi::AllGather( &localNorm, 1, &localNorms[0], 1, rowComm );
    R norm = blas::Nrm2( gridWidth, &localNorms[0], 1 );

    if( norm == 0 )
    {
        if( gridCol == rowAlignment )
            chi.SetLocal(0,0,-chi.GetLocal(0,0));
#ifndef RELEASE
        PopCallStack();
#endif
        return R(2);
    }

    R alpha;
    if( gridCol == rowAlignment )
        alpha = chi.GetLocal(0,0);
    mpi::Broadcast( &alpha, 1, rowAlignment, rowComm );

    R beta;
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

        localNorm = Nrm2( x.LockedMatrix() );
        mpi::AllGather( &localNorm, 1, &localNorms[0], 1, rowComm );
        norm = blas::Nrm2( gridWidth, &localNorms[0], 1 );
        if( alpha <= 0 )
            beta = lapack::SafeNorm( alpha, norm );
        else
            beta = -lapack::SafeNorm( alpha, norm );
    }

    R tau = (beta-alpha)/beta;
    Scale( one/(alpha-beta), x );

    for( int j=0; j<count; ++j )
        beta *= safeInv;
    if( gridCol == rowAlignment )
        chi.SetLocal(0,0,beta);
        
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}

template<typename R>
inline Complex<R>
Row( DistMatrix<Complex<R> >& chi, DistMatrix<Complex<R> >& x )
{
#ifndef RELEASE
    PushCallStack("reflector::Row");    
    if( chi.Grid() != x.Grid() )
        throw std::logic_error
        ("chi and x must be distributed over the same grid");
    if( chi.Height() != 1 || chi.Width() != 1 )
        throw std::logic_error("chi must be a scalar");
    if( x.Height() != 1 )
        throw std::logic_error("x must be a row vector");
    if( chi.Grid().Row() != chi.ColAlignment() )
        throw std::logic_error("Reflecting with incorrect row of processes");
    if( x.Grid().Row() != x.ColAlignment() )
        throw std::logic_error("Reflecting with incorrect row of processes");
#endif
    typedef Complex<R> C;
    const Grid& grid = x.Grid();
    mpi::Comm rowComm = grid.RowComm();
    const int gridCol = grid.Col();
    const int gridWidth = grid.Width();
    const int rowAlignment = chi.RowAlignment();

    std::vector<R> localNorms(gridWidth);
    R localNorm = Nrm2( x.LockedMatrix() ); 
    mpi::AllGather( &localNorm, 1, &localNorms[0], 1, rowComm );
    R norm = blas::Nrm2( gridWidth, &localNorms[0], 1 );

    C alpha;
    if( gridCol == rowAlignment )
        alpha = chi.GetLocal(0,0);
    mpi::Broadcast( &alpha, 1, rowAlignment, rowComm );

    if( norm == R(0) && alpha.imag == R(0) )
    {
        if( gridCol == rowAlignment )
            chi.SetLocal(0,0,-chi.GetLocal(0,0));
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

        localNorm = Nrm2( x.LockedMatrix() );
        mpi::AllGather( &localNorm, 1, &localNorms[0], 1, rowComm );
        norm = blas::Nrm2( gridWidth, &localNorms[0], 1 );
        if( alpha.real <= 0 )
            beta = lapack::SafeNorm( alpha.real, alpha.imag, norm );
        else
            beta = -lapack::SafeNorm( alpha.real, alpha.imag, norm );
    }

    C tau = C( (beta-alpha.real)/beta, -alpha.imag/beta );
    Scale( one/(alpha-beta), x );

    for( int j=0; j<count; ++j )
        beta *= safeInv;
    if( gridCol == rowAlignment )
        chi.SetLocal(0,0,beta);
        
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}

} // namespace reflector
} // namespace elem

#endif // ifndef LAPACK_REFLECTOR_ROW_HPP
