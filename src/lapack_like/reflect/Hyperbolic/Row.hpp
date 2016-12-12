/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HYPERBOLICREFLECTOR_ROW_HPP
#define EL_HYPERBOLICREFLECTOR_ROW_HPP

namespace El {
namespace hyp_reflector {

template<typename F>
F Row( F& chi, AbstractDistMatrix<F>& x )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
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

    vector<Real> localNorms(rowStride);
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
        x *= 1/kappa;
        Conjugate( x );
        return (delta+alpha*lambda)/(kappa*kappa);
    }
}

template<typename F>
F Row( AbstractDistMatrix<F>& chi, AbstractDistMatrix<F>& x )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( chi.ColRank() != chi.ColAlign() || x.ColRank() != x.ColAlign() )
          LogicError("Reflecting from incorrect process");
    )
    F alpha;
    if( chi.IsLocal(0,0) )
        alpha = chi.GetLocal(0,0);
    mpi::Broadcast( alpha, chi.RowAlign(), chi.RowComm() );

    const F tau = reflector::Row( alpha, x );
    chi.Set( 0, 0, alpha );

    return tau;
}

} // namespace hyp_reflector
} // namespace El

#endif // ifndef EL_HYPERBOLICREFLECTOR_ROW_HPP
