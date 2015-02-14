/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// NOTE: This is adapted from a MATLAB script written by AJ Friend.

namespace El {

namespace svm {

template<typename Real>
Int ADMM
( const Matrix<Real>& G, const Matrix<Real>& q, 
        Real gamma,            Matrix<Real>& w,
  Real rho, Int maxIter, bool inv, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("svm::ADMM"))
    const Int numExamples = G.Height();
    const Int numFeatures = G.Width();
    // A = [repmat(q,1,n).*G,q]
    // TODO: Add repmat support into Elemental? It's in Julia.
    Matrix<Real> A( numExamples, numFeatures+1 );
    auto AL = A( IR(0,numExamples), IR(0,numFeatures) );
    AL = G;
    for( Int j=0; j<numFeatures; ++j )
        for( Int i=0; i<numExamples; ++i )
            AL.Set( i, j, AL.Get(i,j)*q.Get(i,0) );
    for( Int i=0; i<numExamples; ++i )
        A.Set( i, numFeatures, q.Get(i,0) );

    auto hingeProx = [=]( Matrix<Real>& y, Real rho )
                     { HingeLossProx( y, y.Height()*rho ); };
    auto frobProx = 
        [=]( Matrix<Real>& x, Real rho ) 
        { auto xT = x( IR(0,x.Height()-1), IR(0,1) );
          FrobeniusProx( xT, gamma/rho ); };

    Matrix<Real> b;
    Zeros( b, numExamples, 1 );

    return ModelFit
    ( function<void(Matrix<Real>&,Real)>(hingeProx),
      function<void(Matrix<Real>&,Real)>(frobProx),
      A, b, w, rho, maxIter, inv, progress );
}

template<typename Real>
Int ADMM
( const AbstractDistMatrix<Real>& G, const AbstractDistMatrix<Real>& q, 
        Real gamma,                        AbstractDistMatrix<Real>& w, 
  Real rho, Int maxIter, bool inv, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("svm::ADMM"))
    const Int numExamples = G.Height();
    const Int numFeatures = G.Width();
    // A = [repmat(q,1,n).*G,q]
    // TODO: Add repmat support into Elemental? It's in Julia.
    DistMatrix<Real> A( numExamples, numFeatures+1, G.Grid() );
    auto AL = A( IR(0,numExamples), IR(0,numFeatures) );
    AL = G;
    DistMatrix<Real,MC,STAR> q_MC_STAR(G.Grid());
    q_MC_STAR.AlignWith( A );
    q_MC_STAR = q;
    for( Int jLoc=0; jLoc<AL.LocalWidth(); ++jLoc )
        for( Int iLoc=0; iLoc<AL.LocalHeight(); ++iLoc )
            AL.SetLocal
            ( iLoc, jLoc, AL.GetLocal(iLoc,jLoc)*q_MC_STAR.GetLocal(iLoc,0) );
    if( A.IsLocalCol(numFeatures) )
    {
        const Int jLoc = A.LocalCol(numFeatures);
        for( Int iLoc=0; iLoc<A.LocalHeight(); ++iLoc )
            A.SetLocal( iLoc, jLoc, q_MC_STAR.GetLocal(iLoc,0) );
    }

    auto hingeProx = [=]( DistMatrix<Real>& y, Real rho )
                     { HingeLossProx( y, y.Height()*rho ); };
    auto frobProx =
        [=]( DistMatrix<Real>& x, Real rho )
        { auto xT = x( IR(0,x.Height()-1), IR(0,1) );
          FrobeniusProx( xT, gamma/rho ); };

    DistMatrix<Real> b(G.Grid());
    Zeros( b, numExamples, 1 );

    return ModelFit
    ( function<void(DistMatrix<Real>&,Real)>(hingeProx),
      function<void(DistMatrix<Real>&,Real)>(frobProx),
      A, b, w, rho, maxIter, inv, progress );
}

#define PROTO(Real) \
  template Int ADMM \
  ( const Matrix<Real>& G, const Matrix<Real>& q, \
          Real gamma,            Matrix<Real>& w, \
    Real rho, Int maxIter, bool inv, bool progress ); \
  template Int ADMM \
  ( const AbstractDistMatrix<Real>& G, const AbstractDistMatrix<Real>& q, \
          Real gamma,                        AbstractDistMatrix<Real>& w, \
    Real rho, Int maxIter, bool inv, bool progress );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace svm

} // namepace elem
