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

    Matrix<Real> d; 
    SOCDets( x, d, orders, firstInds );

    auto Ry = y;
    SOCReflect( Ry, orders, firstInds );
    
    SOCDots( x, y, z, orders, firstInds ); 
    Scale( Real(2), z ); 

    Matrix<Real> detRy;
    Hadamard( d, Ry, detRy );
   
    Axpy( Real(-1), detRy, z );
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
    auto& x = *xPtr;
    auto& y = *yPtr;
    auto& z = *zPtr;

    auto ordersPtr = ReadProxy<Int,VC,STAR>(&ordersPre,ctrl); 
    auto& orders = *ordersPtr;

    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& firstInds = *firstIndsPtr;

    DistMatrix<Real,VC,STAR> d(x.Grid()); 
    SOCDets( x, d, orders, firstInds, cutoff );

    auto Ry = y;
    SOCReflect( Ry, orders, firstInds );

    SOCDots( x, y, z, orders, firstInds, cutoff );
    Scale( Real(2), z );

    DistMatrix<Real,VC,STAR> detRy(x.Grid());
    Hadamard( d, Ry, detRy );

    Axpy( Real(-1), detRy, z );
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

    DistMultiVec<Real> d(x.Comm()); 
    SOCDets( x, d, orders, firstInds, cutoff );

    auto Ry = y;
    SOCReflect( Ry, orders, firstInds );

    SOCDots( x, y, z, orders, firstInds, cutoff );
    Scale( Real(2), z );

    DistMultiVec<Real> detRy(x.Comm());
    Hadamard( d, Ry, detRy );

    Axpy( Real(-1), detRy, z );
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
    const DistMultiVec<Int>& firstInds, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
