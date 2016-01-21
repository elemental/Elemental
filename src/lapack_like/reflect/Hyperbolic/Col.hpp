/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HYPERBOLICREFLECTOR_COL_HPP
#define EL_HYPERBOLICREFLECTOR_COL_HPP

namespace El {
namespace hyp_reflector {

// Please see the comments in the sequential implementation for details.

template<typename F> 
F Col( F& chi, ElementalMatrix<F>& x )
{
    DEBUG_ONLY(
        CSE cse("hyp_reflector::Col");
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

    vector<Real> localNorms(colStride);
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
        x *= 1/kappa;
        return (delta+alpha*lambda)/(kappa*kappa);
    }
}

template<typename F> 
F Col( ElementalMatrix<F>& chi, ElementalMatrix<F>& x )
{
    DEBUG_ONLY(
        CSE cse("hyp_reflector::Col");
        if( chi.RowRank() != chi.RowAlign() || x.RowRank() != x.RowAlign() )
            LogicError("Reflecting from incorrect process");
    )
    F alpha;
    if( chi.IsLocal(0,0) )
        alpha = chi.GetLocal(0,0);
    mpi::Broadcast( alpha, chi.ColAlign(), chi.ColComm() );

    const F tau = reflector::Col( alpha, x );
    chi.Set( 0, 0, alpha );

    return tau;
}

} // namespace hyp_reflector
} // namespace El

#endif // ifndef EL_HYPERBOLICREFLECTOR_COL_HPP
