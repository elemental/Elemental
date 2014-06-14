/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include EL_IDENTITY_INC
#include EL_ONES_INC
#include EL_ZEROS_INC

// NOTE: This is adapted from a MATLAB script written by AJ Friend.

namespace El {

template<typename Real>
Int SVM
( const Matrix<Real>& G, const Matrix<Real>& q, Real gamma, Matrix<Real>& w,
  Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, 
  bool inv, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("SVM"))
    const Int numExamples = G.Height();
    const Int numFeatures = G.Width();
    // A = [repmat(q,1,n).*G,q]
    // TODO: Add repmat support into Elemental? It's in Julia.
    Matrix<Real> A( numExamples, numFeatures+1 );
    auto AL = View( A, 0, 0, numExamples, numFeatures );
    AL = G;
    for( Int j=0; j<numFeatures; ++j )
        for( Int i=0; i<numExamples; ++i )
            AL.Set( i, j, AL.Get(i,j)*q.Get(i,0) );
    for( Int i=0; i<numExamples; ++i )
        A.Set( i, numFeatures, q.Get(i,0) );
    return SVM
    ( A, gamma, w, rho, alpha, maxIter, absTol, relTol, inv, progress );
}

template<typename Real>
Int SVM
( const Matrix<Real>& A, Real gamma, Matrix<Real>& w,
  Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, 
  bool inv, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("SVM"))
    const Int m = A.Height();
    const Int n = A.Width();

    Matrix<Real> P;
    if( m >= n )
    {
        Identity( P, n, n );        
        Herk( LOWER, ADJOINT, Real(1), A, Real(1), P );
    }
    else
    {
        Identity( P, m, m );
        Herk( LOWER, NORMAL, Real(1), A, Real(1), P );
    }
    if( inv )
        HPDInverse( LOWER, P );
    else
        Cholesky( LOWER, P ); 

    // Start the SVM
    Int numIter=0;
    Matrix<Real> s, v1, v2, w1, w2, uw, uv;
    Ones( v2, m, 1 );
    Zeros( w2, n, 1 );
    Zeros( uv, m, 1 );
    Zeros( uw, n, 1 );
    while( numIter < maxIter )
    {
        // Project onto A w1 = v1 
        v1 = v2;
        w1 = w2;
        Axpy( Real(-1), uv, v1 );
        Axpy( Real(-1), uw, w1 );
        Gemv( ADJOINT, Real(1), A, v1, Real(1), w1 );
        if( inv )
        {
            s = w1;
            Hemv( LOWER, Real(1), P, s, Real(0), w1 );
        }
        else
        {
            Trsv( LOWER, NORMAL, NON_UNIT, P, w1 );
            Trsv( LOWER, ADJOINT, NON_UNIT, P, w1 );
        }
        Gemv( NORMAL, Real(1), A, w1, Real(0), v1 );

        // Two-norm prox
        w2 = w1;
        Axpy( Real(1), uw, w2 );
        auto w2T = View( w2, 0, 0, n-1, 1 );
        const Real w2TNorm = FrobeniusNorm( w2T );
        if( w2TNorm > gamma/rho )
            Scale( 1-gamma/(rho*w2TNorm), w2T );
        else
            Zero( w2T );
  
        // Hinge-loss prox
        v2 = v1;
        Axpy( Real(1), uv, v2 );
        EntrywiseMap
        ( v2, [=]( Real alpha ) 
              { if(      alpha <= 1-1/(m*rho) ) return alpha + 1/(m*rho);
                else if( alpha <= Real(1)     ) return Real(1);
                else                            return alpha; } );

        // Update dual variables
        Axpy( Real(1), v1, uv );
        Axpy( Real(1), w1, uw );
        Axpy( Real(-1), v2, uv );
        Axpy( Real(-1), w2, uw );
    }
    if( maxIter == numIter )
        std::cout << "SVM failed to converge" << std::endl;
    w = w2;
    return numIter;
}

template<typename Real>
Int SVM
( const DistMatrix<Real>& G, const DistMatrix<Real>& q, Real gamma, 
  DistMatrix<Real>& w,
  Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, 
  bool inv, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("SVM"))
    const Int numExamples = G.Height();
    const Int numFeatures = G.Width();
    // A = [repmat(q,1,n).*G,q]
    // TODO: Add repmat support into Elemental? It's in Julia.
    DistMatrix<Real> A( numExamples, numFeatures+1, G.Grid() );
    auto AL = View( A, 0, 0, numExamples, numFeatures );
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
    return SVM
    ( A, gamma, w, rho, alpha, maxIter, absTol, relTol, inv, progress );
}

template<typename Real>
Int SVM
( const DistMatrix<Real>& A, Real gamma, DistMatrix<Real>& w,
  Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, 
  bool inv, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("SVM"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();

    DistMatrix<Real> P(g);
    if( m >= n )
    {
        Identity( P, n, n );        
        Herk( LOWER, ADJOINT, Real(1), A, Real(1), P );
    }
    else
    {
        Identity( P, m, m );
        Herk( LOWER, NORMAL, Real(1), A, Real(1), P );
    }
    if( inv )
        HPDInverse( LOWER, P );
    else
        Cholesky( LOWER, P ); 

    // Start the SVM
    Int numIter=0;
    DistMatrix<Real> s(g), v1(g), v2(g), w1(g), w2(g), uw(g), uv(g);
    Ones( v2, m, 1 );
    Zeros( w2, n, 1 );
    Zeros( uv, m, 1 );
    Zeros( uw, n, 1 );
    while( numIter < maxIter )
    {
        // Project onto A w1 = v1 
        v1 = v2;
        w1 = w2;
        Axpy( Real(-1), uv, v1 );
        Axpy( Real(-1), uw, w1 );
        Gemv( ADJOINT, Real(1), A, v1, Real(1), w1 );
        if( inv )
        {
            s = w1;
            Hemv( LOWER, Real(1), P, s, Real(0), w1 );
        }
        else
        {
            Trsv( LOWER, NORMAL, NON_UNIT, P, w1 );
            Trsv( LOWER, ADJOINT, NON_UNIT, P, w1 );
        }
        Gemv( NORMAL, Real(1), A, w1, Real(0), v1 );

        // Two-norm prox
        w2 = w1;
        Axpy( Real(1), uw, w2 );
        auto w2T = View( w2, 0, 0, n-1, 1 );
        const Real w2TNorm = FrobeniusNorm( w2T );
        if( w2TNorm > gamma/rho )
            Scale( 1-gamma/(rho*w2TNorm), w2T );
        else
            Zero( w2T );
  
        // Hinge-loss prox
        v2 = v1;
        Axpy( Real(1), uv, v2 );
        EntrywiseMap
        ( v2, [=]( Real alpha ) 
              { if(      alpha <= 1-1/(m*rho) ) return alpha + 1/(m*rho);
                else if( alpha <= Real(1)     ) return Real(1);
                else                            return alpha; } );

        // Update dual variables
        Axpy( Real(1), v1, uv );
        Axpy( Real(1), w1, uw );
        Axpy( Real(-1), v2, uv );
        Axpy( Real(-1), w2, uw );
    }
    if( maxIter == numIter )
        std::cout << "SVM failed to converge" << std::endl;
    w = w2;
    return numIter;
}

#define PROTO(Real) \
  template Int SVM \
  ( const Matrix<Real>& G, const Matrix<Real>& q, Real gamma, \
    Matrix<Real>& w, \
    Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, \
    bool inv, bool progress ); \
  template Int SVM \
  ( const Matrix<Real>& A, Real gamma, \
    Matrix<Real>& w, \
    Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, \
    bool inv, bool progress ); \
  template Int SVM \
  ( const DistMatrix<Real>& G, const DistMatrix<Real>& q, Real gamma, \
    DistMatrix<Real>& w, \
    Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, \
    bool inv, bool progress ); \
  template Int SVM \
  ( const DistMatrix<Real>& A, Real gamma, \
    DistMatrix<Real>& w, \
    Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, \
    bool inv, bool progress );

PROTO(float)
PROTO(double)

} // namepace elem
