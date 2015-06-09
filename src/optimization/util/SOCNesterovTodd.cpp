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

template<typename Real>
void SOCNesterovTodd
( const Matrix<Real>& s, 
  const Matrix<Real>& z,
        Matrix<Real>& w,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds )
{
    DEBUG_ONLY(CSE cse("SOCNesterovTodd"))

    Matrix<Real> zRoot;
    SOCSquareRoot( z, zRoot, orders, firstInds );

    // a := Q_{sqrt(z)}(s)
    // -------------------
    Matrix<Real> a;
    SOCApplyQuadratic( zRoot, s, a, orders, firstInds );

    // a := sqrt(inv(a)) = inv(sqrt((Q_{sqrt(z)}(s))))
    // -----------------------------------------------
    Matrix<Real> b; 
    SOCInverse( a, b, orders, firstInds );
    SOCSquareRoot( b, a, orders, firstInds );

    // w := Q_{sqrt(z)}(a)
    // -------------------
    SOCApplyQuadratic( zRoot, a, w, orders, firstInds );
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

    auto sPtr = ReadProxy<Real,VC,STAR>(&sPre,ctrl); 
    auto zPtr = ReadProxy<Real,VC,STAR>(&zPre,ctrl);
    auto wPtr = WriteProxy<Real,VC,STAR>(&wPre,ctrl);
    auto ordersPtr = ReadProxy<Int,VC,STAR>(&ordersPre,ctrl); 
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& s = *sPtr;
    auto& z = *zPtr;
    auto& w = *wPtr;
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    DistMatrix<Real,VC,STAR> zRoot(z.Grid());
    SOCSquareRoot( z, zRoot, orders, firstInds, cutoff );

    // a := Q_{sqrt(z)}(s)
    // -------------------
    DistMatrix<Real,VC,STAR> a(s.Grid());
    SOCApplyQuadratic( zRoot, s, a, orders, firstInds, cutoff );

    // a := sqrt(inv(a)) = inv(sqrt((Q_{sqrt(z)}(s))))
    // -----------------------------------------------
    DistMatrix<Real,VC,STAR> b(s.Grid()); 
    SOCInverse( a, b, orders, firstInds, cutoff );
    SOCSquareRoot( b, a, orders, firstInds, cutoff );

    // w := Q_{sqrt(z)}(a)
    // -------------------
    SOCApplyQuadratic( zRoot, a, w, orders, firstInds, cutoff );
}

template<typename Real>
void SOCNesterovTodd
( const DistMultiVec<Real>& s, 
  const DistMultiVec<Real>& z,
        DistMultiVec<Real>& w,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff )
{
    DEBUG_ONLY(CSE cse("SOCNesterovTodd"))

    DistMultiVec<Real> zRoot(z.Comm());
    SOCSquareRoot( z, zRoot, orders, firstInds, cutoff );

    // a := Q_{sqrt(z)}(s)
    // -------------------
    DistMultiVec<Real> a(s.Comm());
    SOCApplyQuadratic( zRoot, s, a, orders, firstInds, cutoff );

    // a := sqrt(inv(a)) = inv(sqrt((Q_{sqrt(z)}(s))))
    // -----------------------------------------------
    DistMultiVec<Real> b(s.Comm()); 
    SOCInverse( a, b, orders, firstInds, cutoff );
    SOCSquareRoot( b, a, orders, firstInds, cutoff );

    // w := Q_{sqrt(z)}(a)
    // -------------------
    SOCApplyQuadratic( zRoot, a, w, orders, firstInds, cutoff );
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
