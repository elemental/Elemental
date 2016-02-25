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

// TODO: Lower-level access

// NOTE: It would be possible to recompute the Nesterov-Todd scaling only for
//       the (almost certainly small) number of modified subcones, but this 
//       would require a non-trivial amount of additional code for a relatively
//       rare and comparatively inexpensive operation.

template<typename Real,typename>
void PushPairInto
(       Matrix<Real>& s, 
        Matrix<Real>& z,
  const Matrix<Real>& w,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds,
  Real wMaxNormLimit )
{
    DEBUG_ONLY(CSE cse("soc::PushPairInto"))

    Matrix<Real> sLower, zLower;
    soc::LowerNorms( s, sLower, orders, firstInds );
    soc::LowerNorms( z, zLower, orders, firstInds );

    const Int height = s.Height();
    for( Int i=0; i<height; ++i )
    {
        const Real w0 = w.Get(i,0);
        const Int firstInd = firstInds.Get(i,0);
        if( i == firstInd && w0 > wMaxNormLimit )
        {
            // TODO: Switch to a non-adhoc modification     
            s.Update( i, 0, Real(1)/wMaxNormLimit );
            z.Update( i, 0, Real(1)/wMaxNormLimit );
        }
    }
}

template<typename Real,typename>
void PushPairInto
(       ElementalMatrix<Real>& sPre, 
        ElementalMatrix<Real>& zPre,
  const ElementalMatrix<Real>& wPre,
  const ElementalMatrix<Int>& ordersPre, 
  const ElementalMatrix<Int>& firstIndsPre,
  Real wMaxNormLimit, Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::PushPairInto"))
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
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = s.GlobalRow(iLoc);
        const Real w0 = w.GetLocal(iLoc,0);
        if( i == firstInds.GetLocal(iLoc,0) && w0 > wMaxNormLimit )
        {
            // TODO: Switch to a non-adhoc modification     
            s.UpdateLocal( iLoc, 0, Real(1)/wMaxNormLimit );
            z.UpdateLocal( iLoc, 0, Real(1)/wMaxNormLimit );
        }
    }
}

template<typename Real,typename>
void PushPairInto
(       DistMultiVec<Real>& s, 
        DistMultiVec<Real>& z,
  const DistMultiVec<Real>& w,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, 
  Real wMaxNormLimit, Int cutoff )
{
    DEBUG_ONLY(CSE cse("soc::PushPairInto"))

    DistMultiVec<Real> sLower(s.Comm()), zLower(z.Comm());
    soc::LowerNorms( s, sLower, orders, firstInds, cutoff );
    soc::LowerNorms( z, zLower, orders, firstInds, cutoff );

    const int localHeight = s.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = s.GlobalRow(iLoc);
        const Real w0 = w.GetLocal(iLoc,0);
        if( i == firstInds.GetLocal(iLoc,0) && w0 > wMaxNormLimit )
        {
            // TODO: Switch to a non-adhoc modification     
            s.UpdateLocal( iLoc, 0, Real(1)/wMaxNormLimit );
            z.UpdateLocal( iLoc, 0, Real(1)/wMaxNormLimit );
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
  (       ElementalMatrix<Real>& s, \
          ElementalMatrix<Real>& z, \
    const ElementalMatrix<Real>& w, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, \
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
#include "El/macros/Instantiate.h"

} // namespace soc
} // namespace El
