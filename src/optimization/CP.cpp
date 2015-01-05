/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// A Chebyshev point (CP) minimizes the supremum of A x - b, i.e.,
//
//   min || A x - b ||_oo.
//
// Real instances of the problem are expressable as a Linear Program [1] via 
//
//   min t
//   s.t. -t <= A x - b <= t,
//
// which, in affine standard form, becomes
//
//   min [1; 0]^T [t; x]
//   s.t. |  A  -1 | | x | <= |  b |
//        | -A  -1 | | t |    | -b |
//
// NOTE: There is likely an appropriate citation, but the derivation is 
//       trivial. If one is found, it will be added.

namespace El {

template<typename Real>
void CP
( const Matrix<Real>& A, const Matrix<Real>& b, 
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("CP"))
    const Int m = A.Height();
    const Int n = A.Width();
    Matrix<Real> c, xHat, AHat, bHat, G, h;

    // c := [zeros(n,1);1]
    // ===================
    Zeros( c, n+1, 1 );
    c.Set( n, 0, Real(1) );

    // \hat A := zeros(0,n+1)
    // ====================== 
    Zeros( AHat, 0, n+1 );

    // \hat b := zeros(0,1)
    // ====================
    Zeros( bHat, 0, 1 );

    // G := |  A -ones(m,1) |
    //      | -A -ones(m,1) |
    // ===========================
    Zeros( G, 2*m, n+1 );
    auto GTL = G( IR(0,m),   IR(0,n)   );
    auto GBL = G( IR(m,2*m), IR(0,n)   );
    auto gTR = G( IR(0,m),   IR(n,n+1) );
    auto gBR = G( IR(m,2*m), IR(n,n+1) );
    GTL = A;
    Axpy( Real(-1), GTL, GBL );
    Fill( gTR, Real(-1) );
    Fill( gBR, Real(-1) );

    // h := |  b |
    //      | -b |
    // ===========
    Zeros( h, 2*m, 1 );
    auto hT = h( IR(0,m),   IR(0,1) );
    auto hB = h( IR(m,2*m), IR(0,1) );
    hT = b;
    hB = b; Scale( Real(-1), hB );

    // Solve the affine LP
    // ===================
    Matrix<Real> y, z, s;
    LP( AHat, G, bHat, c, h, xHat, y, z, s, ctrl );

    // Extract x from [x;t]
    // ====================
    x = xHat( IR(0,n), IR(0,1) );
}

template<typename Real>
void CP
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, 
        AbstractDistMatrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("CP"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& g = A.Grid();
    DistMatrix<Real> c(g), xHat(g), AHat(g), bHat(g), G(g), h(g);

    // c := [zeros(n,1);1]
    // ===================
    Zeros( c, n+1, 1 );
    c.Set( n, 0, Real(1) );

    // \hat A := zeros(0,n+1)
    // ====================== 
    Zeros( AHat, 0, n+1 );

    // \hat b := zeros(0,1)
    // ====================
    Zeros( bHat, 0, 1 );

    // G := |  A -ones(m,1) |
    //      | -A -ones(m,1) |
    // ===========================
    Zeros( G, 2*m, n+1 );
    auto GTL = G( IR(0,m),   IR(0,n)   );
    auto GBL = G( IR(m,2*m), IR(0,n)   );
    auto gTR = G( IR(0,m),   IR(n,n+1) );
    auto gBR = G( IR(m,2*m), IR(n,n+1) );
    GTL = A;
    Axpy( Real(-1), GTL, GBL );
    Fill( gTR, Real(-1) );
    Fill( gBR, Real(-1) );

    // h := |  b |
    //      | -b |
    // ===========
    Zeros( h, 2*m, 1 );
    auto hT = h( IR(0,m),   IR(0,1) );
    auto hB = h( IR(m,2*m), IR(0,1) );
    hT = b;
    hB = b; Scale( Real(-1), hB );

    // Solve the affine LP
    // ===================
    DistMatrix<Real> y(g), z(g), s(g);
    LP( AHat, G, bHat, c, h, xHat, y, z, s, ctrl );

    // Extract x from [x;t]
    // ====================
    Copy( xHat( IR(0,n), IR(0,1) ), x );
}

template<typename Real>
void CP
( const SparseMatrix<Real>& A, const Matrix<Real>& b, 
        Matrix<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("CP"))
    LogicError("This routine is not yet written");
}

template<typename Real>
void CP
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, 
        DistMultiVec<Real>& x,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("CP"))
    LogicError("This routine is not yet written");
}

#define PROTO(Real) \
  template void CP \
  ( const Matrix<Real>& A, const Matrix<Real>& b, \
          Matrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void CP \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, \
          AbstractDistMatrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void CP \
  ( const SparseMatrix<Real>& A, const Matrix<Real>& b, \
          Matrix<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void CP \
  ( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, \
          DistMultiVec<Real>& x, \
    const lp::affine::Ctrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namepace elem
