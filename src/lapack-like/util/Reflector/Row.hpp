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
#ifndef EL_REFLECTOR_ROW_HPP
#define EL_REFLECTOR_ROW_HPP

namespace El {
namespace reflector {

template<typename F,Dist U,Dist V>
F Row( DistMatrix<F,U,V>& chi, DistMatrix<F,U,V>& x )
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

    F tau = (beta-Conj(alpha)) / beta;
    Scale( Real(1)/(alpha-beta), x );

    // Undo the scaling
    for( Int j=0; j<count; ++j )
        beta *= safeInv;

    if( rowRank == rowAlign )
        chi.SetLocal(0,0,beta);

    // This is to make this a reflector meant to be applied from the right;
    // there is no need to conjugate chi, as it is real
    Conjugate( x );
        
    return tau;
}

template<typename F,Dist U,Dist V>
F Row( F& chi, DistMatrix<F,U,V>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("reflector::Row");    
        if( x.Height() != 1 )
            LogicError("x must be a row vector");
        if( x.ColRank() != x.ColAlign() )
            LogicError("Reflecting from incorrect process");
    )
    typedef Base<F> Real;
    mpi::Comm rowComm = x.RowComm();
    const Int rowStride = x.RowStride();

    std::vector<Real> localNorms(rowStride);
    Real localNorm = Nrm2( x.LockedMatrix() ); 
    mpi::AllGather( &localNorm, 1, localNorms.data(), 1, rowComm );
    Real norm = blas::Nrm2( rowStride, localNorms.data(), 1 );

    F alpha = chi;
    if( norm == Real(0) && ImagPart(alpha) == Real(0) )
    {
        chi = -chi;
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

    F tau = (beta-Conj(alpha)) / beta;
    Scale( Real(1)/(alpha-beta), x );

    // Undo the scaling
    for( Int j=0; j<count; ++j )
        beta *= safeInv;

    chi = beta;

    // This is to make this a reflector meant to be applied from the right;
    // there is no need to conjugate chi, as it is real
    Conjugate( x );
        
    return tau;
}

} // namespace reflector
} // namespace El

#endif // ifndef EL_REFLECTOR_ROW_HPP
