/*
   Copyright (c) 2009-2016, Jack Poulson
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
  const ModelFitCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("svm::ADMM"))
    const Int numExamples = G.Height();
    const Int numFeatures = G.Width();
    // A = [repmat(q,1,n).*G,q]
    // TODO: Add repmat support into Elemental? It's in Julia.
    Matrix<Real> A( numExamples, numFeatures+1 );
    auto AL = A( ALL, IR(0,numFeatures) );
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
        { auto xT = x( IR(0,x.Height()-1), ALL );
          FrobeniusProx( xT, gamma/rho ); };

    Matrix<Real> b;
    Zeros( b, numExamples, 1 );

    return ModelFit
    ( function<void(Matrix<Real>&,Real)>(hingeProx),
      function<void(Matrix<Real>&,Real)>(frobProx),
      A, b, w, ctrl );
}

template<typename Real>
Int ADMM
( const ElementalMatrix<Real>& G, const ElementalMatrix<Real>& q, 
        Real gamma,                        ElementalMatrix<Real>& w, 
  const ModelFitCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("svm::ADMM"))
    const Int numExamples = G.Height();
    const Int numFeatures = G.Width();
    // A = [repmat(q,1,n).*G,q]
    // TODO: Add repmat support into Elemental? It's in Julia.
    DistMatrix<Real> A( numExamples, numFeatures+1, G.Grid() );
    auto AL = A( ALL, IR(0,numFeatures) );
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
        { auto xT = x( IR(0,x.Height()-1), ALL );
          FrobeniusProx( xT, gamma/rho ); };

    DistMatrix<Real> b(G.Grid());
    Zeros( b, numExamples, 1 );

    return ModelFit
    ( function<void(DistMatrix<Real>&,Real)>(hingeProx),
      function<void(DistMatrix<Real>&,Real)>(frobProx),
      A, b, w, ctrl );
}

} // namespace svm
} // namespace El
