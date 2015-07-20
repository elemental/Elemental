/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Q_x y = (2 x x^T - det(x) R) y = 2 (x^T y) x - det(x) (R y)

template<typename Real>
void SOCApplyQuadratic
( const Matrix<Real>& x, 
  const Matrix<Real>& y,
        Matrix<Real>& z,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("SOCApplyQuadratic"))

    // detRy := det(x) R y
    Matrix<Real> d;
    SOCDets( x, d, orders, firstInds );
    ConeBroadcast( d, orders, firstInds );
    auto Ry = y;
    SOCReflect( Ry, orders, firstInds );
    Matrix<Real> detRy;
    Hadamard( d, Ry, detRy );

    // z := 2 (x^T y) x
    Matrix<Real> xTy;
    SOCDots( x, y, xTy, orders, firstInds );
    ConeBroadcast( xTy, orders, firstInds );
    Hadamard( xTy, x, z );
    z *= 2;

    // z := z - detRy
    z -= detRy;
}

template<typename Real>
void SOCApplyQuadratic
( const AbstractDistMatrix<Real>& xPre, 
  const AbstractDistMatrix<Real>& yPre,
        AbstractDistMatrix<Real>& zPre,
  const AbstractDistMatrix<Int>& ordersPre, 
  const AbstractDistMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCApplyQuadratic"))
    AssertSameGrids( xPre, yPre, zPre, ordersPre, firstIndsPre );

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    auto xPtr = ReadProxy<Real,VC,STAR>(&xPre,ctrl); 
    auto yPtr = ReadProxy<Real,VC,STAR>(&yPre,ctrl);
    auto zPtr = WriteProxy<Real,VC,STAR>(&zPre,ctrl);
    auto ordersPtr = ReadProxy<Int,VC,STAR>(&ordersPre,ctrl); 
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& x = *xPtr;
    auto& y = *yPtr;
    auto& z = *zPtr;
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    // detRy := det(x) R y
    DistMatrix<Real,VC,STAR> d(x.Grid());
    SOCDets( x, d, orders, firstInds, cutoff );
    ConeBroadcast( d, orders, firstInds, cutoff );
    auto Ry = y;
    SOCReflect( Ry, orders, firstInds );
    DistMatrix<Real,VC,STAR> detRy(x.Grid());
    Hadamard( d, Ry, detRy );

    // z := 2 (x^T y) x
    DistMatrix<Real,VC,STAR> xTy(x.Grid());
    SOCDots( x, y, xTy, orders, firstInds, cutoff );
    ConeBroadcast( xTy, orders, firstInds, cutoff );
    Hadamard( xTy, x, z );
    z *= 2;

    // z := z - detRy
    z -= detRy;
}

template<typename Real>
void SOCApplyQuadratic
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCApplyQuadratic"))

    // detRy := det(x) R y
    DistMultiVec<Real> d(x.Comm()); 
    SOCDets( x, d, orders, firstInds, cutoff );
    ConeBroadcast( d, orders, firstInds, cutoff );
    auto Ry = y;
    SOCReflect( Ry, orders, firstInds );
    DistMultiVec<Real> detRy(x.Comm());
    Hadamard( d, Ry, detRy );

    // z := 2 (x^T y) x
    DistMultiVec<Real> xTy(x.Comm());
    SOCDots( x, y, xTy, orders, firstInds, cutoff );
    ConeBroadcast( xTy, orders, firstInds, cutoff );
    Hadamard( xTy, x, z );
    z *= 2;

    // z := z - detRy
    z -= detRy;
}

template<typename Real>
void SOCApplyQuadratic
( const Matrix<Real>& x,
        Matrix<Real>& y,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("SOCApplyQuadratic"))
    // TODO?: Optimize
    Matrix<Real> z; 
    SOCApplyQuadratic( x, y, z, orders, firstInds );
    y = z;
}

template<typename Real>
void SOCApplyQuadratic
( const AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& y,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCApplyQuadratic"))
    // TODO?: Optimize
    DistMatrix<Real,VC,STAR> z(x.Grid()); 
    SOCApplyQuadratic( x, y, z, orders, firstInds, cutoff );
    y = z;
}

template<typename Real>
void SOCApplyQuadratic
( const DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCApplyQuadratic"))
    // TODO?: Optimize 
    DistMultiVec<Real> z(x.Comm());
    SOCApplyQuadratic( x, y, z, orders, firstInds, cutoff );
    y = z;
}

#define PROTO(Real) \
  template void SOCApplyQuadratic \
  ( const Matrix<Real>& x, \
    const Matrix<Real>& y, \
          Matrix<Real>& z, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void SOCApplyQuadratic \
  ( const AbstractDistMatrix<Real>& x, \
    const AbstractDistMatrix<Real>& y, \
          AbstractDistMatrix<Real>& z, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, Int cutoff ); \
  template void SOCApplyQuadratic \
  ( const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff ); \
  template void SOCApplyQuadratic \
  ( const Matrix<Real>& x, \
          Matrix<Real>& y, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void SOCApplyQuadratic \
  ( const AbstractDistMatrix<Real>& x, \
          AbstractDistMatrix<Real>& y, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, Int cutoff ); \
  template void SOCApplyQuadratic \
  ( const DistMultiVec<Real>& x, \
          DistMultiVec<Real>& y, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
