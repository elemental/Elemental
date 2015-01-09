/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// The soft-margin Support Vector Machine (SVM) [1] seeks the solution of the
// Quadratic Program
//
//   min_{w,beta,z} (1/2) || w ||_2^2 + lambda 1^T z
//
//   s.t. |-diag(d) A, -d, -I | | w    | <= | -1 |
//        |      0,     0, -I | | beta |    |  0 |
//                              | z    |
//
// [1] Corinna Cortes and Vladimir Vapnik,
//     "Support-Vector Networks",
//     Journal of Machine Learning, Vol. 20, No. 3, 1995.
//

namespace El {

template<typename Real>
void SVM
( const Matrix<Real>& A, const Matrix<Real>& d, 
        Real lambda,           Matrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("SVM"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Range<Int> wInd(0,n), betaInd(n,n+1), zInd(n+1,n+m+1);

    Matrix<Real> Q, c, AHat, b, G, h;

    // Q := | I 0 0 |
    //      | 0 0 0 |
    //      | 0 0 0 |
    // ==============
    Zeros( Q, n+m+1, n+m+1 );
    auto Qww = Q( wInd, wInd );
    FillDiagonal( Qww, Real(1) );

    // c := [0;0;lambda]
    // =================
    Zeros( c, n+m+1, 1 );
    auto cz = c( zInd, IR(0,1) );
    Fill( cz, lambda );

    // AHat = []
    // =========
    Zeros( AHat, 0, n+m+1 );

    // b = []
    // ======
    Zeros( b, 0, 1 );

    // G := |-diag(d) A, -d, -I|
    //      |      0,     0, -I|
    // =========================
    Zeros( G, 2*m, n+m+1 );
    auto G0w    = G( IR(0,m),   wInd    );
    auto G0beta = G( IR(0,m),   betaInd );
    auto G0z    = G( IR(0,m),   zInd    );
    auto G1z    = G( IR(m,2*m), zInd    );
    G0w = A; Scale( Real(-1), G0w ); DiagonalScale( LEFT, NORMAL, d, G0w );
    G0beta = d; Scale( Real(-1), G0beta );
    FillDiagonal( G0z, Real(-1) );
    FillDiagonal( G1z, Real(-1) );

    // h := [-ones(m,1); zeros(m,1)]
    // =============================
    Zeros( h, 2*m, 1 );
    auto h0 = h( IR(0,m), IR(0,1) );
    Fill( h0, Real(-1) );
    
    // Solve the affine QP
    // ===================
    Matrix<Real> y, z, s;
    QP( Q, AHat, G, b, c, h, x, y, z, s, ctrl );
}

template<typename Real>
void SVM
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& d, 
        Real lambda,                       AbstractDistMatrix<Real>& x, 
  const qp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("svm::ADMM"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    const Range<Int> wInd(0,n), betaInd(n,n+1), zInd(n+1,n+m+1);

    DistMatrix<Real> Q(g), c(g), AHat(g), b(g), G(g), h(g);

    // Q := | I 0 0 |
    //      | 0 0 0 |
    //      | 0 0 0 |
    // ==============
    Zeros( Q, n+m+1, n+m+1 );
    auto Qww = Q( wInd, wInd );
    FillDiagonal( Qww, Real(1) );

    // c := [0;0;lambda]
    // ================
    Zeros( c, n+m+1, 1 );
    auto cz = c( zInd, IR(0,1) );
    Fill( cz, lambda );

    // AHat = []
    // =========
    Zeros( AHat, 0, n+m+1 );

    // b = []
    // ======
    Zeros( b, 0, 1 );

    // G := |-diag(d) A, -d, -I|
    //      |      0,     0, -I|
    // =========================
    Zeros( G, 2*m, n+m+1 );
    auto G0w    = G( IR(0,m),   wInd    );
    auto G0beta = G( IR(0,m),   betaInd );
    auto G0z    = G( IR(0,m),   zInd    );
    auto G1z    = G( IR(m,2*m), zInd    );
    G0w = A; Scale( Real(-1), G0w ); DiagonalScale( LEFT, NORMAL, d, G0w );
    G0beta = d; Scale( Real(-1), G0beta );
    FillDiagonal( G0z, Real(-1) );
    FillDiagonal( G1z, Real(-1) );

    // h := [-ones(m,1); zeros(m,1)]
    // =============================
    Zeros( h, 2*m, 1 );
    auto h0 = h( IR(0,m), IR(0,1) );
    Fill( h0, Real(-1) );
    
    // Solve the affine QP
    // ===================
    DistMatrix<Real> y(g), z(g), s(g);
    QP( Q, AHat, G, b, c, h, x, y, z, s, ctrl );
}

template<typename Real>
void SVM
( const SparseMatrix<Real>& A, const Matrix<Real>& d, 
        Real lambda,                 Matrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("SVM"))
    LogicError("This routine is not yet written");
}

template<typename Real>
void SVM
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& d, 
        Real lambda,                     DistMultiVec<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("SVM"))
    LogicError("This routine is not yet written");
}

#define PROTO(Real) \
  template void SVM \
  ( const Matrix<Real>& A, const Matrix<Real>& d, \
          Real lambda,           Matrix<Real>& x, \
    const qp::affine::Ctrl<Real>& ctrl ); \
  template void SVM \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& d, \
          Real lambda,                       AbstractDistMatrix<Real>& x, \
    const qp::affine::Ctrl<Real>& ctrl ); \
  template void SVM \
  ( const SparseMatrix<Real>& A, const Matrix<Real>& d, \
          Real lambda,                 Matrix<Real>& x, \
    const qp::affine::Ctrl<Real>& ctrl ); \
  template void SVM \
  ( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& d, \
          Real lambda,                     DistMultiVec<Real>& x, \
    const qp::affine::Ctrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namepace elem
