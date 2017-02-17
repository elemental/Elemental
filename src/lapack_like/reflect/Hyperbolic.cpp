/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./Hyperbolic/Col.hpp"
#include "./Hyperbolic/Row.hpp"

namespace El {

//
// (I - 1/tau w w^H Sigma) x = -lambda e_0
// 
// where Sigma is diag(+1,-1,...,-1)

template<typename F>
F LeftHyperbolicReflector( Matrix<F>& chi, Matrix<F>& x )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( chi.Height() != 1 || chi.Width() != 1 )
          LogicError("chi must be a scalar");
    )

    F alpha = chi(0);
    const F tau = LeftHyperbolicReflector( alpha, x );
    chi(0) = alpha;

    return tau;
}

template<typename F>
F LeftHyperbolicReflector( F& chi, Matrix<F>& x )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( x.Height() != 1 && x.Width() != 1 )
          LogicError("x must be a vector");
      if( ImagPart(chi) != Base<F>(0) )
          LogicError("chi is assumed to be real");
    )

    // Compute lambda = sgn(chi) sqrt([chi;x]^H Sigma [chi;x])
    //                = sgn(chi) sqrt(chi^2 - x^H x)    
    typedef Base<F> Real;
    const Real alpha = RealPart(chi);
    const Real xNrm = Nrm2( x ); 
    const Real delta = alpha*alpha - xNrm*xNrm;
    if( delta < Real(0) )
        LogicError("Attempted to square-root a negative number");
    const Real lambda = ( alpha>=0 ? Sqrt(delta) : -Sqrt(delta) );
    chi = -lambda;

    // Implicitly define 
    //     w := [chi;x] + lambda e_0, and 
    //     kappa = chi + lambda,
    // so that
    //     tau := (w^H Sigma w) / 2 = (delta + lambda^2 + 2 chi lambda) / 2
    //          = delta + chi lambda
    // then normalize w so that its first entry is one
    // TODO: Introduce a threshold instead of the approach from 
    //       van de Geijn and van Zee's "High-performance up-and-downdating
    //       via Householder-like transformations"
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
F LeftHyperbolicReflector
( AbstractDistMatrix<F>& chi, AbstractDistMatrix<F>& x )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( chi, x );
      if( chi.Height() != 1 || chi.Width() != 1 )
          LogicError("chi must be a scalar");
      if( x.Width() != 1 )
          LogicError("x must be a column vector");
      if( chi.Root() != x.Root() )
          LogicError("Roots must be the same");
    )
    F tau;
    if( x.CrossRank() == x.Root() )
    {
        if( x.RowRank() == x.RowAlign() )
            tau = hyp_reflector::Col( chi, x );
        mpi::Broadcast( tau, x.RowAlign(), x.RowComm() );
    }
    mpi::Broadcast( tau, x.Root(), x.CrossComm() );
    return tau;
}

template<typename F>
F LeftHyperbolicReflector( F& chi, AbstractDistMatrix<F>& x )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( x.Width() != 1 )
          LogicError("x must be a column vector");
    )
    F tau;
    if( x.CrossRank() == x.Root() )
    {
        if( x.RowRank() == x.RowAlign() )
            tau = hyp_reflector::Col( chi, x );
        mpi::Broadcast( tau, x.RowAlign(), x.RowComm() );
    }
    mpi::Broadcast( tau, x.Root(), x.CrossComm() );
    return tau;
}

//
// x (I - 1/tau Sigma w^H w) = -lambda e_0^T
// 
// where Sigma is diag(+1,-1,...,-1)

template<typename F>
F RightHyperbolicReflector( Matrix<F>& chi, Matrix<F>& x )
{
    EL_DEBUG_CSE
    const F tau = LeftHyperbolicReflector( chi, x );
    Conjugate( x );
    return tau;
}

template<typename F>
F RightHyperbolicReflector( F& chi, Matrix<F>& x )
{
    EL_DEBUG_CSE
    const F tau = LeftHyperbolicReflector( chi, x );
    Conjugate( x );
    return tau;
}

template<typename F>
F RightHyperbolicReflector
( AbstractDistMatrix<F>& chi, AbstractDistMatrix<F>& x )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( chi, x );
      if( chi.Height() != 1 || chi.Width() != 1 )
          LogicError("chi must be a scalar");
      if( x.Height() != 1 )
          LogicError("x must be a row vector");
      if( chi.Root() != x.Root() )
          LogicError("Roots must be the same");
    )
    F tau;
    if( x.CrossRank() == x.Root() )
    {
        if( x.ColRank() == x.ColAlign() )
            tau = hyp_reflector::Row( chi, x );
        mpi::Broadcast( tau, x.ColAlign(), x.ColComm() );
    }
    mpi::Broadcast( tau, x.Root(), x.CrossComm() );
    return tau;
}

template<typename F>
F RightHyperbolicReflector( F& chi, AbstractDistMatrix<F>& x )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( x.Height() != 1 )
          LogicError("x must be a row vector");
    )
    F tau;
    if( x.CrossRank() == x.Root() )
    {
        if( x.ColRank() == x.ColAlign() )
            tau = hyp_reflector::Row( chi, x );
        mpi::Broadcast( tau, x.ColAlign(), x.ColComm() );
    }
    mpi::Broadcast( tau, x.Root(), x.CrossComm() );
    return tau;
}

#define PROTO(F) \
  template F LeftHyperbolicReflector( F& chi, Matrix<F>& x ); \
  template F LeftHyperbolicReflector( F& chi, AbstractDistMatrix<F>& x ); \
  template F LeftHyperbolicReflector( Matrix<F>& chi, Matrix<F>& x ); \
  template F LeftHyperbolicReflector \
  ( AbstractDistMatrix<F>& chi, AbstractDistMatrix<F>& x ); \
  template F RightHyperbolicReflector( F& chi, Matrix<F>& x ); \
  template F RightHyperbolicReflector( F& chi, AbstractDistMatrix<F>& x ); \
  template F RightHyperbolicReflector( Matrix<F>& chi, Matrix<F>& x ); \
  template F RightHyperbolicReflector \
  ( AbstractDistMatrix<F>& chi, AbstractDistMatrix<F>& x ); \
  template F hyp_reflector::Col( F& chi, AbstractDistMatrix<F>& x ); \
  template F hyp_reflector::Col \
  ( AbstractDistMatrix<F>& chi, AbstractDistMatrix<F>& x ); \
  template F hyp_reflector::Row( F& chi, AbstractDistMatrix<F>& x ); \
  template F hyp_reflector::Row \
  ( AbstractDistMatrix<F>& chi, AbstractDistMatrix<F>& x );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
