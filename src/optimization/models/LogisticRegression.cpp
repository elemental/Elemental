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

template<typename Real>
Int LogisticRegression
( const Matrix<Real>& G, const Matrix<Real>& q, Matrix<Real>& w,
  Real gamma, Regularization penalty, 
  const ModelFitCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("LogisticRegression"))
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

    auto logisticProx = [&]( Matrix<Real>& y, Real rho )
                        { LogisticProx( y, y.Height()*rho ); };
    auto logisticFunc = function<void(Matrix<Real>&,Real)>(logisticProx);

    function<void(Matrix<Real>&,Real)> proxFunc;
    if( penalty == NO_PENALTY )
    {
        auto noProx = [&]( Matrix<Real>& x, Real rho ) { };
        proxFunc = function<void(Matrix<Real>&,Real)>(noProx);
    }
    else if( penalty == L1_PENALTY )
    {
        auto oneProx = 
            [&]( Matrix<Real>& x, Real rho ) 
            { auto xT = x( IR(0,x.Height()-1), ALL );
              SoftThreshold( xT, gamma/rho ); };
        proxFunc = function<void(Matrix<Real>&,Real)>(oneProx);
    }
    else if( penalty == L2_PENALTY )
    {
        auto frobProx = 
            [&]( Matrix<Real>& x, Real rho ) 
            { auto xT = x( IR(0,x.Height()-1), ALL );
              FrobeniusProx( xT, gamma/rho ); };
        proxFunc = function<void(Matrix<Real>&,Real)>(frobProx);
    }

    Matrix<Real> b;
    Zeros( b, numExamples, 1 );

    return ModelFit( logisticFunc, proxFunc, A, b, w, ctrl );
}

template<typename Real>
Int LogisticRegression
( const ElementalMatrix<Real>& G, const ElementalMatrix<Real>& q, 
        ElementalMatrix<Real>& w, 
  Real gamma, Regularization penalty, 
  const ModelFitCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("LogisticRegression"))
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

    auto logisticProx = [&]( DistMatrix<Real>& y, Real rho )
                        { LogisticProx( y, y.Height()*rho ); };
    auto logisticFunc = function<void(DistMatrix<Real>&,Real)>(logisticProx);
    
    function<void(DistMatrix<Real>&,Real)> proxFunc;
    if( penalty == NO_PENALTY )
    {
        auto noProx = [&]( DistMatrix<Real>& x, Real rho ) { };
        proxFunc = function<void(DistMatrix<Real>&,Real)>(noProx);
    }
    else if( penalty == L1_PENALTY )
    {    
        auto oneProx =
            [&]( DistMatrix<Real>& x, Real rho )
            { auto xT = x( IR(0,x.Height()-1), ALL );
              SoftThreshold( xT, gamma/rho ); };
        proxFunc = function<void(DistMatrix<Real>&,Real)>(oneProx);
    }
    else if( penalty == L2_PENALTY )
    {    
        auto frobProx =
            [&]( DistMatrix<Real>& x, Real rho )
            { auto xT = x( IR(0,x.Height()-1), ALL );
              FrobeniusProx( xT, gamma/rho ); };
        proxFunc = function<void(DistMatrix<Real>&,Real)>(frobProx);
    }

    DistMatrix<Real> b(G.Grid());
    Zeros( b, numExamples, 1 );

    return ModelFit( logisticFunc, proxFunc, A, b, w, ctrl );
}

#define PROTO(Real) \
  template Int LogisticRegression \
  ( const Matrix<Real>& G, const Matrix<Real>& q, Matrix<Real>& w, \
    Real gamma, Regularization penalty, \
    const ModelFitCtrl<Real>& ctrl ); \
  template Int LogisticRegression \
  ( const ElementalMatrix<Real>& G, const ElementalMatrix<Real>& q, \
          ElementalMatrix<Real>& w, \
    Real gamma, Regularization penalty, \
    const ModelFitCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
