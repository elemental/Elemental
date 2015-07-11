/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "../../affine/IPM.hpp"

namespace El {
namespace socp {
namespace direct {

// The following solves a pair of SOC programs in "direct" conic form:
//
//   min c^T x
//   s.t. A x = b, x in K,
//
//   max -b^T y
//   s.t. A^T y - z + c = 0, z in K,
//
// as opposed to the more general "affine" conic form:
//
//   min c^T x
//   s.t. A x = b, G x + s = h, s in K,
//
//   max -b^T y - h^T z
//   s.t. A^T y + G^T z + c = 0, z in K
//

template<typename Real>
void Mehrotra
( const Matrix<Real>& A, 
  const Matrix<Real>& b, 
  const Matrix<Real>& c,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& x, 
        Matrix<Real>& y, 
        Matrix<Real>& z, 
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("socp::direct::Mehrotra"))    
    const Int n = c.Height();

    Matrix<Real> G;
    Identity( G, n, n );
    G *= -1;

    Matrix<Real> h;
    Zeros( h, n, 1 );

    MehrotraCtrl<Real> affineCtrl;
    affineCtrl.primalInit = false;
    affineCtrl.dualInit = false;
    affineCtrl.minTol = ctrl.minTol;
    affineCtrl.targetTol = ctrl.targetTol;
    affineCtrl.maxIts = ctrl.maxIts;
    affineCtrl.maxStepRatio = ctrl.maxStepRatio;
    affineCtrl.qsdCtrl = ctrl.qsdCtrl;
    affineCtrl.outerEquil = ctrl.outerEquil;
    affineCtrl.innerEquil = ctrl.innerEquil;
    affineCtrl.basisSize = ctrl.basisSize;
    affineCtrl.print = ctrl.print;
    affineCtrl.time = ctrl.time;

    Matrix<Real> s;
    socp::affine::Mehrotra(A,G,b,c,h,orders,firstInds,x,y,z,s,affineCtrl); 
}

template<typename Real>
void Mehrotra
( const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& b, 
  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
        AbstractDistMatrix<Real>& x, 
        AbstractDistMatrix<Real>& y, 
        AbstractDistMatrix<Real>& z, 
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("socp::direct::Mehrotra"))    
    const Int n = c.Height();
    const Grid& grid = c.Grid();

    DistMatrix<Real> G(grid); 
    Identity( G, n, n );
    G *= -1;

    DistMatrix<Real> h(grid);
    Zeros( h, n, 1 );

    MehrotraCtrl<Real> affineCtrl;
    affineCtrl.primalInit = false;
    affineCtrl.dualInit = false;
    affineCtrl.minTol = ctrl.minTol;
    affineCtrl.targetTol = ctrl.targetTol;
    affineCtrl.maxIts = ctrl.maxIts;
    affineCtrl.maxStepRatio = ctrl.maxStepRatio;
    affineCtrl.qsdCtrl = ctrl.qsdCtrl;
    affineCtrl.outerEquil = ctrl.outerEquil;
    affineCtrl.innerEquil = ctrl.innerEquil;
    affineCtrl.basisSize = ctrl.basisSize;
    affineCtrl.print = ctrl.print;
    affineCtrl.time = ctrl.time;

    DistMatrix<Real> s(grid);
    socp::affine::Mehrotra(A,G,b,c,h,orders,firstInds,x,y,z,s,affineCtrl);
}

template<typename Real>
void Mehrotra
( const SparseMatrix<Real>& A, 
  const Matrix<Real>& b, 
  const Matrix<Real>& c,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("socp::direct::Mehrotra"))    
    const Int n = c.Height();

    SparseMatrix<Real> G;
    Identity( G, n, n );
    G *= -1;

    Matrix<Real> h;
    Zeros( h, n, 1 );

    MehrotraCtrl<Real> affineCtrl;
    affineCtrl.primalInit = false;
    affineCtrl.dualInit = false;
    affineCtrl.minTol = ctrl.minTol;
    affineCtrl.targetTol = ctrl.targetTol;
    affineCtrl.maxIts = ctrl.maxIts;
    affineCtrl.maxStepRatio = ctrl.maxStepRatio;
    affineCtrl.qsdCtrl = ctrl.qsdCtrl;
    affineCtrl.outerEquil = ctrl.outerEquil;
    affineCtrl.innerEquil = ctrl.innerEquil;
    affineCtrl.basisSize = ctrl.basisSize;
    affineCtrl.print = ctrl.print;
    affineCtrl.time = ctrl.time;

    Matrix<Real> s;
    socp::affine::Mehrotra(A,G,b,c,h,orders,firstInds,x,y,z,s,affineCtrl);
}

template<typename Real>
void Mehrotra
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("socp::direct::Mehrotra"))    
    const Int n = c.Height();
    mpi::Comm comm = c.Comm();

    DistSparseMatrix<Real> G(comm);
    Identity( G, n, n );
    G *= -1;

    DistMultiVec<Real> h(comm);
    Zeros( h, n, 1 );

    MehrotraCtrl<Real> affineCtrl;
    affineCtrl.primalInit = false;
    affineCtrl.dualInit = false;
    affineCtrl.minTol = ctrl.minTol;
    affineCtrl.targetTol = ctrl.targetTol;
    affineCtrl.maxIts = ctrl.maxIts;
    affineCtrl.maxStepRatio = ctrl.maxStepRatio;
    affineCtrl.qsdCtrl = ctrl.qsdCtrl;
    affineCtrl.outerEquil = ctrl.outerEquil;
    affineCtrl.innerEquil = ctrl.innerEquil;
    affineCtrl.basisSize = ctrl.basisSize;
    affineCtrl.print = ctrl.print;
    affineCtrl.time = ctrl.time;

    DistMultiVec<Real> s(comm);
    socp::affine::Mehrotra(A,G,b,c,h,orders,firstInds,x,y,z,s,affineCtrl);
}

#define PROTO(Real) \
  template void Mehrotra \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b, \
    const AbstractDistMatrix<Real>& c, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, \
          AbstractDistMatrix<Real>& x, \
          AbstractDistMatrix<Real>& y, \
          AbstractDistMatrix<Real>& z, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, \
    const DistMultiVec<Real>& c, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
          DistMultiVec<Real>& x, \
          DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
    const MehrotraCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace direct
} // namespace socp
} // namespace El
