/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HYPERBOLICREFLECTOR_COL_HPP
#define ELEM_HYPERBOLICREFLECTOR_COL_HPP

#include ELEM_NRM2_INC
#include ELEM_SCALE_INC
#include ELEM_ZERO_INC

namespace elem {
namespace hyp_reflector {

// Please see the comments in the sequential implementation for details.

template<typename F,Dist U,Dist V> 
inline F
Col( DistMatrix<F,U,V>& chi, DistMatrix<F,U,V>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("hyp_reflector::Col");
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

    Real alpha;
    if( colRank == colAlign )
    {
        if( ImagPart(chi.GetLocal(0,0)) != Real(0) )
            LogicError("chi is assumed to be real");
        alpha = chi.GetLocalRealPart(0,0);
    }
    mpi::Broadcast( alpha, colAlign, colComm );
    const Real delta = alpha*alpha - norm*norm;
    if( delta < Real(0) )
        LogicError("Attempted to square-root a negative number");
    const Real lambda = ( alpha>=0 ? Sqrt(delta) : -Sqrt(delta) );
    if( colRank == colAlign ) 
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
        return (delta+alpha*lambda)/(kappa*kappa);
    }
}

template<typename F,Dist U,Dist V> 
inline F
Col( F& chi, DistMatrix<F,U,V>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("hyp_reflector::Col");
        if( x.Width() != 1 )
            LogicError("x must be a column vector");
        if( x.RowRank() != x.RowAlign() )
            LogicError("Reflecting from incorrect process");
        if( ImagPart(chi) !=Base<F>(0) )
            LogicError("chi is assumed to be real");
    )
    typedef Base<F> Real;
    mpi::Comm colComm = x.ColComm();
    const Int colStride = x.ColStride();

    std::vector<Real> localNorms(colStride);
    Real localNorm = Nrm2( x.LockedMatrix() ); 
    mpi::AllGather( &localNorm, 1, localNorms.data(), 1, colComm );
    Real norm = blas::Nrm2( colStride, localNorms.data(), 1 );

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
        return (delta+alpha*lambda)/(kappa*kappa);
    }
}

} // namespace hyp_reflector
} // namespace elem

#endif // ifndef ELEM_HYPERBOLICREFLECTOR_COL_HPP
