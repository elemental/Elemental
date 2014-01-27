/*
   Copyright (c) 1992-2008 The University of Tennessee. 
   All rights reserved.

   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is loosely based upon the LAPACK routines dlarfg.f and zlarfg.f.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_REFLECTOR_ROW_HPP
#define ELEM_REFLECTOR_ROW_HPP

#include ELEM_NRM2_INC
#include ELEM_SCALE_INC

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

template<typename F,Dist U,Dist V>
inline F
Row( DistMatrix<F,U,V>& chi, DistMatrix<F,U,V>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("reflector::Row");    
        if( chi.Grid() != x.Grid() )
            LogicError("chi and x must be distributed over the same grid");
        if( chi.Height() != 1 || chi.Width() != 1 )
            LogicError("chi must be a scalar");
        if( x.Height() != 1 )
            LogicError("x must be a row vector");
        if( chi.ColRank() != chi.ColAlign() || x.ColRank() != x.ColAlign() )
            LogicError("Reflecting from incorrect process");
    )
    typedef Base<F> Real;
    mpi::Comm rowComm = x.RowComm();
    const Int rowRank = x.RowRank();
    const Int rowStride = x.RowStride();
    const Int rowAlign = chi.RowAlign();

    std::vector<Real> localNorms(rowStride);
    Real localNorm = Nrm2( x.LockedMatrix() ); 
    mpi::AllGather( &localNorm, 1, localNorms.data(), 1, rowComm );
    Real norm = blas::Nrm2( rowStride, localNorms.data(), 1 );

    F alpha;
    if( rowRank == rowAlign )
        alpha = chi.GetLocal(0,0);
    mpi::Broadcast( alpha, rowAlign, rowComm );

    if( norm == Real(0) && ImagPart(alpha) == Real(0) )
    {
        if( rowRank == rowAlign )
            chi.SetLocal(0,0,-chi.GetLocal(0,0));
        return F(2);
    }

    Real beta;
    if( RealPart(alpha) <= 0 )
        beta = lapack::SafeNorm( alpha, norm );
    else
        beta = -lapack::SafeNorm( alpha, norm );

    // Rescale if the vector is too small
    const Real safeMin = lapack::MachineSafeMin<Real>();
    const Real epsilon = lapack::MachineEpsilon<Real>();
    const Real safeInv = safeMin/epsilon;
    Int count = 0;
    if( Abs(beta) < safeInv )
    {
        Real invOfSafeInv = Real(1)/safeInv;
        do
        {
            ++count;
            Scale( invOfSafeInv, x );
            alpha *= invOfSafeInv;
            beta *= invOfSafeInv;
        } while( Abs(beta) < safeInv );

        localNorm = Nrm2( x.LockedMatrix() );
        mpi::AllGather( &localNorm, 1, localNorms.data(), 1, rowComm );
        norm = blas::Nrm2( rowStride, localNorms.data(), 1 );
        if( RealPart(alpha) <= 0 )
            beta = lapack::SafeNorm( alpha, norm );
        else
            beta = -lapack::SafeNorm( alpha, norm );
    }

    F tau = (beta-alpha) / beta;
    Scale( Real(1)/(alpha-beta), x );

    // Undo the scaling
    for( Int j=0; j<count; ++j )
        beta *= safeInv;

    if( rowRank == rowAlign )
        chi.SetLocal(0,0,beta);
        
    return tau;
}

} // namespace reflector
} // namespace elem

#endif // ifndef ELEM_REFLECTOR_ROW_HPP
