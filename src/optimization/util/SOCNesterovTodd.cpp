/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// NOTE: Many of the determinant calculations happening below are redundant.

// Find the Nesterov-Todd scaling point w such that 
//
//   Q_w z = s,
//
// or, equivalently, for each subcone,
//
//   Hess(-(1/2) ln det(w)) s = z,
//
// where -(1/2) ln det(w) is the barrier function for the Second-Order Cone.

template<typename Real>
void SOCNesterovTodd
( const Matrix<Real>& s, 
  const Matrix<Real>& z,
        Matrix<Real>& w,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("SOCNesterovTodd"))

    Matrix<Promote<Real>> sProm, zProm;
    Copy( s, sProm );
    Copy( z, zProm );

    Matrix<Promote<Real>> sRoot;
    SOCSquareRoot( sProm, sRoot, orders, firstInds );

    // a := Q_{sqrt(s)}(z)
    // -------------------
    Matrix<Promote<Real>> a;
    SOCApplyQuadratic( sRoot, zProm, a, orders, firstInds );

    // a := sqrt(inv(a)) = inv(sqrt((Q_{sqrt(s)}(z))))
    // -----------------------------------------------
    Matrix<Promote<Real>> b; 
    SOCInverse( a, b, orders, firstInds );
    SOCSquareRoot( b, a, orders, firstInds );

    // w := Q_{sqrt(s)}(a)
    // -------------------
    Matrix<Promote<Real>> wProm;
    SOCApplyQuadratic( sRoot, a, wProm, orders, firstInds );
    Copy( wProm, w );
}

template<typename Real>
void SOCNesterovTodd
( const AbstractDistMatrix<Real>& sPre, 
  const AbstractDistMatrix<Real>& zPre,
        AbstractDistMatrix<Real>& wPre,
  const AbstractDistMatrix<Int>& ordersPre, 
  const AbstractDistMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCNesterovTodd"))
    AssertSameGrids( sPre, zPre, wPre, ordersPre, firstIndsPre );

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    auto sPtr = ReadProxy<Promote<Real>,VC,STAR>(&sPre,ctrl); 
    auto zPtr = ReadProxy<Promote<Real>,VC,STAR>(&zPre,ctrl);
    auto wPtr = WriteProxy<Promote<Real>,VC,STAR>(&wPre,ctrl);
    auto ordersPtr = ReadProxy<Int,VC,STAR>(&ordersPre,ctrl); 
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& s = *sPtr;
    auto& z = *zPtr;
    auto& w = *wPtr;
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    DistMatrix<Promote<Real>,VC,STAR> sRoot(s.Grid());
    SOCSquareRoot( s, sRoot, orders, firstInds, cutoff );

    // a := Q_{sqrt(s)}(z)
    // -------------------
    DistMatrix<Promote<Real>,VC,STAR> a(z.Grid());
    SOCApplyQuadratic( sRoot, z, a, orders, firstInds, cutoff );

    // a := sqrt(inv(a)) = inv(sqrt((Q_{sqrt(s)}(z))))
    // -----------------------------------------------
    DistMatrix<Promote<Real>,VC,STAR> b(z.Grid()); 
    SOCInverse( a, b, orders, firstInds, cutoff );
    SOCSquareRoot( b, a, orders, firstInds, cutoff );

    // w := Q_{sqrt(s)}(a)
    // -------------------
    SOCApplyQuadratic( sRoot, a, w, orders, firstInds, cutoff );
}

// TODO: Extend promotions to other cases
template<typename Real>
void SOCNesterovTodd
( const DistMultiVec<Real>& s, 
  const DistMultiVec<Real>& z,
        DistMultiVec<Real>& w,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCNesterovTodd"))
    mpi::Comm comm = s.Comm();

    DistMultiVec<Promote<Real>> sProm(comm), zProm(comm);
    Copy( s, sProm );
    Copy( z, zProm );

    DistMultiVec<Promote<Real>> sRoot(comm);
    SOCSquareRoot( sProm, sRoot, orders, firstInds, cutoff );

    // a := Q_{sqrt(s)}(z)
    // -------------------
    DistMultiVec<Promote<Real>> a(comm);
    SOCApplyQuadratic( sRoot, zProm, a, orders, firstInds, cutoff );

    // a := sqrt(inv(a)) = inv(sqrt((Q_{sqrt(s)}(z))))
    // -----------------------------------------------
    DistMultiVec<Promote<Real>> b(comm); 
    SOCInverse( a, b, orders, firstInds, cutoff );
    SOCSquareRoot( b, a, orders, firstInds, cutoff );

    // w := Q_{sqrt(s)}(a)
    // -------------------
    DistMultiVec<Promote<Real>> wProm(comm);
    SOCApplyQuadratic( sRoot, a, wProm, orders, firstInds, cutoff );
    Copy( wProm, w );
}

#define PROTO(Real) \
  template void SOCNesterovTodd \
  ( const Matrix<Real>& s, \
    const Matrix<Real>& z, \
          Matrix<Real>& w, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void SOCNesterovTodd \
  ( const AbstractDistMatrix<Real>& s, \
    const AbstractDistMatrix<Real>& z, \
          AbstractDistMatrix<Real>& w, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, Int cutoff ); \
  template void SOCNesterovTodd \
  ( const DistMultiVec<Real>& s, \
    const DistMultiVec<Real>& z, \
          DistMultiVec<Real>& w, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
