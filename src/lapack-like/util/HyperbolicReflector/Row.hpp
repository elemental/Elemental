/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_HYPERBOLICREFLECTOR_ROW_HPP
#define EL_HYPERBOLICREFLECTOR_ROW_HPP

namespace El {
namespace hyp_reflector {

template<typename F,Dist U,Dist V>
F Row( DistMatrix<F,U,V>& chi, DistMatrix<F,U,V>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("hyp_reflector::Row");
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

    Real alpha;
    if( rowRank == rowAlign )
    {
        if( ImagPart(chi.GetLocal(0,0)) != Real(0) )
            LogicError("chi is assumed to be real");
        alpha = chi.GetLocalRealPart(0,0);
    }
    mpi::Broadcast( alpha, rowAlign, rowComm );
    const Real delta = alpha*alpha - norm*norm;
    if( delta < Real(0) )
        LogicError("Attempted to square-root a negative number");
    const Real lambda = ( alpha>=0 ? Sqrt(delta) : -Sqrt(delta) );
    if( rowRank == rowAlign )
        chi.SetLocal(0,0,-lambda);

    const Real kappa = alpha + lambda;
    if( kappa == Real(0) )
    {
        Zero( x );
        return Real(1);
    }
    else
    {
        Scale( Real(1)/kappa, x );
        Conjugate( x );
        return (delta+alpha*lambda)/(kappa*kappa);
    }
}

template<typename F,Dist U,Dist V>
F Row( F& chi, DistMatrix<F,U,V>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("hyp_reflector::Row");
        if( x.Height() != 1 )
            LogicError("x must be a row vector");
        if( x.ColRank() != x.ColAlign() )
            LogicError("Reflecting from incorrect process");
        if( ImagPart(chi) != Base<F>(0) )
            LogicError("chi is assumed to be real");
    )
    typedef Base<F> Real;
    mpi::Comm rowComm = x.RowComm();
    const Int rowStride = x.RowStride();

    std::vector<Real> localNorms(rowStride);
    Real localNorm = Nrm2( x.LockedMatrix() );
    mpi::AllGather( &localNorm, 1, localNorms.data(), 1, rowComm );
    Real norm = blas::Nrm2( rowStride, localNorms.data(), 1 );

    const Real alpha = RealPart(chi);
    const Real delta = alpha*alpha - norm*norm;
    if( delta < Real(0) )
        LogicError("Attempted to square-root a negative number");
    const Real lambda = ( alpha>=0 ? Sqrt(delta) : -Sqrt(delta) );
    chi = -lambda;

    const Real kappa = alpha + lambda;
    if( kappa == Real(0) )
    {
        Zero( x );
        return Real(1);
    }
    else
    {
        Scale( Real(1)/kappa, x );
        Conjugate( x );
        return (delta+alpha*lambda)/(kappa*kappa);
    }
}

} // namespace hyp_reflector
} // namespace El

#endif // ifndef EL_HYPERBOLICREFLECTOR_ROW_HPP
