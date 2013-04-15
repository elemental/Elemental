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
#ifndef LAPACK_REFLECTOR_COL_HPP
#define LAPACK_REFLECTOR_COL_HPP

#include "elemental/blas-like/level1/Nrm2.hpp"
#include "elemental/blas-like/level1/Scale.hpp"

namespace elem {
namespace reflector {

//
// Follows the LAPACK convention of defining tau such that
//
//   H = I - tau [1; v] [1, v'],
//
// but adjoint(H) [chi; x] = [beta; 0]. 
//
// Note that the adjoint of H is applied. In the case of real data,
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
Col( DistMatrix<R>& chi, DistMatrix<R>& x )
{
#ifndef RELEASE
    PushCallStack("reflector::Col");
    if( chi.Grid() != x.Grid() )
        throw std::logic_error
        ("chi and x must be distributed over the same grid");
    if( chi.Height() != 1 || chi.Width() != 1 )
        throw std::logic_error("chi must be a scalar");
    if( x.Width() != 1 )
        throw std::logic_error("x must be a column vector");
    if( chi.Grid().Col() != chi.RowAlignment() )
        throw std::logic_error("Reflecting with incorrect column of processes");
    if( x.Grid().Col() != x.RowAlignment() )
        throw std::logic_error("Reflecting with incorrect column of processes");
#endif
    const Grid& grid = x.Grid();
    mpi::Comm colComm = grid.ColComm();
    const int gridHeight = grid.Height();
    const int gridRow = grid.Row();
    const int colAlignment = chi.ColAlignment();

    std::vector<R> localNorms(gridHeight);
    R localNorm = Nrm2( x.LockedMatrix() ); 
    mpi::AllGather( &localNorm, 1, &localNorms[0], 1, colComm );
    R norm = blas::Nrm2( gridHeight, &localNorms[0], 1 );

    if( norm == 0 )
    {
        if( gridRow == colAlignment )
            chi.SetLocal(0,0,-chi.GetLocal(0,0));
#ifndef RELEASE
        PopCallStack();
#endif
        return R(2);
    }

    R alpha;
    if( gridRow == colAlignment )
        alpha = chi.GetLocal(0,0);
    mpi::Broadcast( &alpha, 1, colAlignment, colComm );

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
        mpi::AllGather( &localNorm, 1, &localNorms[0], 1, colComm );
        norm = blas::Nrm2( gridHeight, &localNorms[0], 1 );
        if( alpha <= 0 )
            beta = lapack::SafeNorm( alpha, norm );
        else
            beta = -lapack::SafeNorm( alpha, norm );
    }

    R tau = (beta-alpha)/beta;
    Scale( one/(alpha-beta), x );

    for( int j=0; j<count; ++j )
        beta *= safeInv;
    if( gridRow == colAlignment )
        chi.SetLocal(0,0,beta);
        
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}

template<typename R> 
inline Complex<R>
Col( DistMatrix<Complex<R> >& chi, DistMatrix<Complex<R> >& x )
{
#ifndef RELEASE
    PushCallStack("reflector::Col");
    if( chi.Grid() != x.Grid() )
        throw std::logic_error
        ("chi and x must be distributed over the same grid");
    if( chi.Height() != 1 || chi.Width() != 1 )
        throw std::logic_error("chi must be a scalar");
    if( x.Width() != 1 )
        throw std::logic_error("x must be a column vector");
    if( chi.Grid().Col() != chi.RowAlignment() )
        throw std::logic_error("Reflecting with incorrect column of processes");
    if( x.Grid().Col() != x.RowAlignment() )
        throw std::logic_error("Reflecting with incorrect column of processes");
#endif
    typedef Complex<R> C;
    const Grid& grid = x.Grid();
    mpi::Comm colComm = grid.ColComm();
    const int gridHeight = grid.Height();
    const int gridRow = grid.Row();
    const int colAlignment = chi.ColAlignment();

    std::vector<R> localNorms(gridHeight);
    R localNorm = Nrm2( x.LockedMatrix() ); 
    mpi::AllGather( &localNorm, 1, &localNorms[0], 1, colComm );
    R norm = blas::Nrm2( gridHeight, &localNorms[0], 1 );

    C alpha;
    if( gridRow == colAlignment )
        alpha = chi.GetLocal(0,0);
    mpi::Broadcast( &alpha, 1, colAlignment, colComm );

    if( norm == R(0) && alpha.imag == R(0) )
    {
        if( gridRow == colAlignment )
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
        mpi::AllGather( &localNorm, 1, &localNorms[0], 1, colComm );
        norm = blas::Nrm2( gridHeight, &localNorms[0], 1 );
        if( alpha.real <= 0 )
            beta = lapack::SafeNorm( alpha.real, alpha.imag, norm );
        else
            beta = -lapack::SafeNorm( alpha.real, alpha.imag, norm );
    }

    C tau = C( (beta-alpha.real)/beta, -alpha.imag/beta );
    Scale( one/(alpha-beta), x );

    for( int j=0; j<count; ++j )
        beta *= safeInv;
    if( gridRow == colAlignment )
        chi.SetLocal(0,0,beta);
        
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}

} // namespace reflector
} // namespace elem

#endif // ifndef LAPACK_REFLECTOR_COL_HPP
