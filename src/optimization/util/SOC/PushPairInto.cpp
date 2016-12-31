/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace soc {

// TODO(poulson): Lower-level access

// NOTE: It would be possible to recompute the Nesterov-Todd scaling only for
//       the (almost certainly small) number of modified subcones, but this
//       would require a non-trivial amount of additional code for a relatively
//       rare and comparatively inexpensive operation.

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void PushPairInto
(       Matrix<Real>& s,
        Matrix<Real>& z,
  const Matrix<Real>& w,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  Real wMaxNormLimit )
{
    EL_DEBUG_CSE

    Matrix<Real> sLower, zLower;
    soc::LowerNorms( s, sLower, orders, firstInds );
    soc::LowerNorms( z, zLower, orders, firstInds );

    const Int height = s.Height();
    for( Int i=0; i<height; ++i )
    {
        if( i == firstInds(i) && w(i) > wMaxNormLimit )
        {
            // TODO(poulson): Switch to a non-adhoc modification
            s(i) += Real(1)/wMaxNormLimit;
            z(i) += Real(1)/wMaxNormLimit;
        }
    }
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void PushPairInto
(       AbstractDistMatrix<Real>& sPre,
        AbstractDistMatrix<Real>& zPre,
  const AbstractDistMatrix<Real>& wPre,
  const AbstractDistMatrix<Int>& ordersPre,
  const AbstractDistMatrix<Int>& firstIndsPre,
  Real wMaxNormLimit, Int cutoff )
{
    EL_DEBUG_CSE
    AssertSameGrids( sPre, zPre, wPre, ordersPre, firstIndsPre );

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    DistMatrixReadWriteProxy<Real,Real,VC,STAR>
      sProx( sPre, ctrl ),
      zProx( zPre, ctrl );
    DistMatrixReadProxy<Real,Real,VC,STAR>
      wProx( wPre, ctrl );
    DistMatrixReadProxy<Int,Int,VC,STAR>
      ordersProx( ordersPre, ctrl ),
      firstIndsProx( firstIndsPre, ctrl );
    auto& s = sProx.Get();
    auto& z = zProx.Get();
    auto& w = wProx.GetLocked();
    auto& orders = ordersProx.GetLocked();
    auto& firstInds = firstIndsProx.GetLocked();

    DistMatrix<Real,VC,STAR> sLower(s.Grid()), zLower(z.Grid());
    soc::LowerNorms( s, sLower, orders, firstInds, cutoff );
    soc::LowerNorms( z, zLower, orders, firstInds, cutoff );

    const Int localHeight = s.LocalHeight();
    auto& sLoc = s.Matrix();
    auto& zLoc = z.Matrix();
    auto& wLoc = w.LockedMatrix();
    auto& firstIndsLoc = firstInds.LockedMatrix();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = s.GlobalRow(iLoc);
        if( i == firstIndsLoc(iLoc) && wLoc(iLoc) > wMaxNormLimit )
        {
            // TODO(poulson): Switch to a non-adhoc modification
            sLoc(iLoc) += Real(1)/wMaxNormLimit;
            zLoc(iLoc) += Real(1)/wMaxNormLimit;
        }
    }
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void PushPairInto
(       DistMultiVec<Real>& s,
        DistMultiVec<Real>& z,
  const DistMultiVec<Real>& w,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Real wMaxNormLimit, Int cutoff )
{
    EL_DEBUG_CSE

    DistMultiVec<Real> sLower(s.Grid()), zLower(z.Grid());
    soc::LowerNorms( s, sLower, orders, firstInds, cutoff );
    soc::LowerNorms( z, zLower, orders, firstInds, cutoff );

    const int localHeight = s.LocalHeight();
    auto& sLoc = s.Matrix();
    auto& zLoc = z.Matrix();
    auto& wLoc = w.LockedMatrix();
    auto& firstIndsLoc = firstInds.LockedMatrix();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = s.GlobalRow(iLoc);
        if( i == firstIndsLoc(iLoc) && wLoc(iLoc) > wMaxNormLimit )
        {
            // TODO(poulson): Switch to a non-adhoc modification
            sLoc(iLoc) += Real(1)/wMaxNormLimit;
            zLoc(iLoc) += Real(1)/wMaxNormLimit;
        }
    }
}

#define PROTO(Real) \
  template void PushPairInto \
  (       Matrix<Real>& s, \
          Matrix<Real>& z, \
    const Matrix<Real>& w, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    Real wMaxNormLimit ); \
  template void PushPairInto \
  (       AbstractDistMatrix<Real>& s, \
          AbstractDistMatrix<Real>& z, \
    const AbstractDistMatrix<Real>& w, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, \
    Real wMaxNormLimit, Int cutoff ); \
  template void PushPairInto \
  (       DistMultiVec<Real>& s, \
          DistMultiVec<Real>& z, \
    const DistMultiVec<Real>& w, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    Real wMaxNormLimit, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace soc
} // namespace El
