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
#ifndef EL_REFLECTOR_COL_HPP
#define EL_REFLECTOR_COL_HPP

namespace El {
namespace reflector {

template<typename F,Dist U,Dist V> 
F Col( DistMatrix<F,U,V>& chi, DistMatrix<F,U,V>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("reflector::Col");
        if( chi.Grid() != x.Grid() )
            LogicError("chi and x must be distributed over the same grid");
        if( chi.Height() != 1 || chi.Width() != 1 )
            LogicError("chi must be a scalar");
        if( x.Width() != 1 )
            LogicError("x must be a column vector");
        if( chi.RowRank() != chi.RowAlign() || x.RowRank() != x.RowAlign() )
            LogicError("Reflecting from incorrect process");
    )
    typedef Base<F> Real;
    mpi::Comm colComm = x.ColComm();
    const Int colRank = x.ColRank();
    const Int colStride = x.ColStride();
    const Int colAlign = chi.ColAlign();

    std::vector<Real> localNorms(colStride);
    Real localNorm = Nrm2( x.LockedMatrix() ); 
    mpi::AllGather( &localNorm, 1, localNorms.data(), 1, colComm );
    Real norm = blas::Nrm2( colStride, localNorms.data(), 1 );

    F alpha;
    if( colRank == colAlign )
        alpha = chi.GetLocal(0,0);
    mpi::Broadcast( alpha, colAlign, colComm );

    if( norm == Real(0) && ImagPart(alpha) == Real(0) )
    {
        if( colRank == colAlign )
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
        mpi::AllGather( &localNorm, 1, localNorms.data(), 1, colComm );
        norm = blas::Nrm2( colStride, localNorms.data(), 1 );
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

    if( colRank == colAlign )
        chi.SetLocal(0,0,beta);
        
    return tau;
}

template<typename F,Dist U,Dist V> 
F Col( F& chi, DistMatrix<F,U,V>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("reflector::Col");
        if( x.Width() != 1 )
            LogicError("x must be a column vector");
        if( x.RowRank() != x.RowAlign() )
            LogicError("Reflecting from incorrect process");
    )
    typedef Base<F> Real;
    mpi::Comm colComm = x.ColComm();
    const Int colStride = x.ColStride();

    std::vector<Real> localNorms(colStride);
    Real localNorm = Nrm2( x.LockedMatrix() ); 
    mpi::AllGather( &localNorm, 1, localNorms.data(), 1, colComm );
    Real norm = blas::Nrm2( colStride, localNorms.data(), 1 );

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
        mpi::AllGather( &localNorm, 1, localNorms.data(), 1, colComm );
        norm = blas::Nrm2( colStride, localNorms.data(), 1 );
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
    return tau;
}

} // namespace reflector
} // namespace El

#endif // ifndef EL_REFLECTOR_COL_HPP
