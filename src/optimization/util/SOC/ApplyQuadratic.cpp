/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace soc {

// Q_x y = (2 x x^T - det(x) R) y = 2 (x^T y) x - det(x) (R y)

template<typename Real,typename>
void ApplyQuadratic
( const Matrix<Real>& x, 
  const Matrix<Real>& y,
        Matrix<Real>& z,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("soc::ApplyQuadratic"))

    // detRy := det(x) R y
    Matrix<Real> d;
    soc::Dets( x, d, orders, firstInds );
    cone::Broadcast( d, orders, firstInds );
    auto Ry = y;
    soc::Reflect( Ry, orders, firstInds );
    Matrix<Real> detRy;
    Hadamard( d, Ry, detRy );

    // z := 2 (x^T y) x
    Matrix<Real> xTy;
    soc::Dots( x, y, xTy, orders, firstInds );
    cone::Broadcast( xTy, orders, firstInds );
    Hadamard( xTy, x, z );
    z *= 2;

    // z := z - detRy
    z -= detRy;
}

template<typename Real,typename>
void ApplyQuadratic
( const ElementalMatrix<Real>& xPre, 
  const ElementalMatrix<Real>& yPre,
        ElementalMatrix<Real>& zPre,
  const ElementalMatrix<Int>& ordersPre, 
  const ElementalMatrix<Int>& firstIndsPre,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::ApplyQuadratic"))
    AssertSameGrids( xPre, yPre, zPre, ordersPre, firstIndsPre );

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    DistMatrixReadProxy<Real,Real,VC,STAR>
      xProx( xPre, ctrl ),
      yProx( yPre, ctrl );
    DistMatrixWriteProxy<Real,Real,VC,STAR>
      zProx( zPre, ctrl );
    DistMatrixReadProxy<Int,Int,VC,STAR>
      ordersProx( ordersPre, ctrl ),
      firstIndsProx( firstIndsPre, ctrl );
    auto& x = xProx.GetLocked();
    auto& y = yProx.GetLocked();
    auto& z = zProx.Get();
    auto& orders = ordersProx.GetLocked();
    auto& firstInds = firstIndsProx.GetLocked(); 

    // detRy := det(x) R y
    DistMatrix<Real,VC,STAR> d(x.Grid());
    soc::Dets( x, d, orders, firstInds, cutoff );
    cone::Broadcast( d, orders, firstInds, cutoff );
    auto Ry = y;
    soc::Reflect( Ry, orders, firstInds );
    DistMatrix<Real,VC,STAR> detRy(x.Grid());
    Hadamard( d, Ry, detRy );

    // z := 2 (x^T y) x
    DistMatrix<Real,VC,STAR> xTy(x.Grid());
    soc::Dots( x, y, xTy, orders, firstInds, cutoff );
    cone::Broadcast( xTy, orders, firstInds, cutoff );
    Hadamard( xTy, x, z );
    z *= 2;

    // z := z - detRy
    z -= detRy;
}

template<typename Real,typename>
void ApplyQuadratic
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::ApplyQuadratic"))

    // detRy := det(x) R y
    DistMultiVec<Real> d(x.Comm()); 
    soc::Dets( x, d, orders, firstInds, cutoff );
    cone::Broadcast( d, orders, firstInds, cutoff );
    auto Ry = y;
    soc::Reflect( Ry, orders, firstInds );
    DistMultiVec<Real> detRy(x.Comm());
    Hadamard( d, Ry, detRy );

    // z := 2 (x^T y) x
    DistMultiVec<Real> xTy(x.Comm());
    soc::Dots( x, y, xTy, orders, firstInds, cutoff );
    cone::Broadcast( xTy, orders, firstInds, cutoff );
    Hadamard( xTy, x, z );
    z *= 2;

    // z := z - detRy
    z -= detRy;
}

template<typename Real,typename>
void ApplyQuadratic
( const Matrix<Real>& x,
        Matrix<Real>& y,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("soc::ApplyQuadratic"))
    // TODO?: Optimize
    Matrix<Real> z; 
    soc::ApplyQuadratic( x, y, z, orders, firstInds );
    y = z;
}

template<typename Real,typename>
void ApplyQuadratic
( const ElementalMatrix<Real>& x,
        ElementalMatrix<Real>& y,
  const ElementalMatrix<Int>& orders,
  const ElementalMatrix<Int>& firstInds,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::ApplyQuadratic"))
    // TODO?: Optimize
    DistMatrix<Real,VC,STAR> z(x.Grid()); 
    soc::ApplyQuadratic( x, y, z, orders, firstInds, cutoff );
    y = z;
}

template<typename Real,typename>
void ApplyQuadratic
( const DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::ApplyQuadratic"))
    // TODO?: Optimize 
    DistMultiVec<Real> z(x.Comm());
    soc::ApplyQuadratic( x, y, z, orders, firstInds, cutoff );
    y = z;
}

#define PROTO(Real) \
  template void ApplyQuadratic \
  ( const Matrix<Real>& x, \
    const Matrix<Real>& y, \
          Matrix<Real>& z, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void ApplyQuadratic \
  ( const ElementalMatrix<Real>& x, \
    const ElementalMatrix<Real>& y, \
          ElementalMatrix<Real>& z, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, \
    Int cutoff ); \
  template void ApplyQuadratic \
  ( const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    Int cutoff ); \
  template void ApplyQuadratic \
  ( const Matrix<Real>& x, \
          Matrix<Real>& y, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds ); \
  template void ApplyQuadratic \
  ( const ElementalMatrix<Real>& x, \
          ElementalMatrix<Real>& y, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, \
    Int cutoff ); \
  template void ApplyQuadratic \
  ( const DistMultiVec<Real>& x, \
          DistMultiVec<Real>& y, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace soc
} // namespace El
