/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include EL_IDENTITY_INC

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
    Matrix<Real> s, x0, x1, x2, y0, y1, y2, ux, uy;

    Zeros( x1, n, 1 );
    Zeros( x2, n, 1 );
    Zeros( ux, n, 1 );

    Zeros( y1, m, 1 );
    Ones( y2, m, 1 );
    Zeros( uy, m, 1 );

    for( ; numIter<maxIter-1; ++numIter )
    //while( numIter < maxIter )
    {
        // Project onto A x1 + b = y1
        x0 = x2; 
        y0 = y2;
        Axpy( Real(-1), ux, x0 );
        Axpy( Real(-1), uy, y0 );
        x1 = x0;
        // TODO: Subtract A' b
        Gemv( ADJOINT, Real(1), A, y0, Real(1), x1 );
        if( inv )
        {
            s = x1;
            Hemv( LOWER, Real(1), P, s, Real(0), x1 );
        }
        else
        {
            Trsv( LOWER, NORMAL, NON_UNIT, P, x1 );
            Trsv( LOWER, ADJOINT, NON_UNIT, P, x1 );
        }
        Gemv( NORMAL, Real(1), A, x1, Real(0), y1 );

        // Hinge-loss prox
        y2 = y1;
        Axpy( Real(1), uy, y2 );
        EntrywiseMap
        ( y2, [=]( Real alpha ) 
              { if( alpha < 1 ) { return Min(alpha+1/(m*rho),Real(1)); } 
                else              return alpha;                          } );

        x2 = x1; 
        Axpy( Real(1), ux, x2 );
        // TODO: Turn this two-norm prox into a subroutine
        auto x2T = View( x2, 0, 0, n-1, 1 );
        const Real x2TNorm = FrobeniusNorm( x2T );
        if( x2TNorm > gamma/rho )
            Scale( 1-gamma/(rho*x2TNorm), x2T );
        else
            Zero( x2T );

        // Update dual variables
        Axpy( Real(1), x1, ux );
        Axpy( Real(1), y1, uy );
        Axpy( Real(-1), x2, ux );
        Axpy( Real(-1), y2, uy );
    }
    if( maxIter == numIter )
        std::cout << "SVM failed to converge" << std::endl;
    w = x2;
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
    DistMatrix<Real> s(g), x0(g), x1(g), x2(g), y0(g), y1(g), y2(g), 
                     ux(g), uy(g);

    Zeros( x1, n, 1 );
    Zeros( x2, n, 1 );
    Zeros( ux, n, 1 );

    Zeros( y1, m, 1 );
    Ones( y2, m, 1 );
    Zeros( uy, m, 1 );

    for( ; numIter<maxIter-1; ++numIter )
    //while( numIter < maxIter )
    {
        // Project onto A x1 + b = y1
        x0 = x2; 
        y0 = y2;
        Axpy( Real(-1), ux, x0 );
        Axpy( Real(-1), uy, y0 );
        x1 = x0;
        // TODO: Subtract A' b
        Gemv( ADJOINT, Real(1), A, y0, Real(1), x1 );
        if( inv )
        {
            s = x1;
            Hemv( LOWER, Real(1), P, s, Real(0), x1 );
        }
        else
        {
            Trsv( LOWER, NORMAL, NON_UNIT, P, x1 );
            Trsv( LOWER, ADJOINT, NON_UNIT, P, x1 );
        }
        Gemv( NORMAL, Real(1), A, x1, Real(0), y1 );

        // Hinge-loss prox
        y2 = y1;
        Axpy( Real(1), uy, y2 );
        EntrywiseMap
        ( y2, [=]( Real alpha ) 
              { if( alpha < 1 ) { return Min(alpha+1/(m*rho),Real(1)); } 
                else              return alpha;                          } );

        x2 = x1; 
        Axpy( Real(1), ux, x2 );
        // TODO: Turn this two-norm prox into a subroutine
        auto x2T = View( x2, 0, 0, n-1, 1 );
        const Real x2TNorm = FrobeniusNorm( x2T );
        if( x2TNorm > gamma/rho )
            Scale( 1-gamma/(rho*x2TNorm), x2T );
        else
            Zero( x2T );

        // Update dual variables
        Axpy( Real(1), x1, ux );
        Axpy( Real(1), y1, uy );
        Axpy( Real(-1), x2, ux );
        Axpy( Real(-1), y2, uy );
    }
    if( maxIter == numIter )
        std::cout << "SVM failed to converge" << std::endl;
    w = x2;
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
