/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// NOTE: It would be possible to recompute the Nesterov-Todd scaling only for
//       the (almost certainly small) number of modified subcones, but this 
//       would require a non-trivial amount of additional code for a relatively
//       rare and comparatively inexpensive operation.

template<typename Real>
void ForcePairIntoSOC
(       Matrix<Real>& s, 
        Matrix<Real>& z,
  const Matrix<Real>& w,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds,
  Real wMaxNormLimit )
{
    DEBUG_ONLY(CSE cse("ForcePairIntoSOC"))

    Matrix<Real> sLower, zLower;
    SOCLowerNorms( s, sLower, orders, firstInds );
    SOCLowerNorms( z, zLower, orders, firstInds );

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

template<typename Real>
void ForcePairIntoSOC
(       AbstractDistMatrix<Real>& sPre, 
        AbstractDistMatrix<Real>& zPre,
  const AbstractDistMatrix<Real>& wPre,
  const AbstractDistMatrix<Int>& ordersPre, 
  const AbstractDistMatrix<Int>& firstIndsPre,
  Real wMaxNormLimit, Int cutoff )
{
    DEBUG_ONLY(CSE cse("ForcePairIntoSOC"))
    AssertSameGrids( sPre, zPre, wPre, ordersPre, firstIndsPre );

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    auto sPtr = ReadWriteProxy<Real,VC,STAR>(&sPre,ctrl); 
    auto zPtr = ReadWriteProxy<Real,VC,STAR>(&zPre,ctrl); 
    auto wPtr = ReadProxy<Real,VC,STAR>(&wPre,ctrl);
    auto ordersPtr = ReadProxy<Int,VC,STAR>(&ordersPre,ctrl); 
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre,ctrl);
    auto& s = *sPtr;
    auto& z = *zPtr;
    auto& w = *wPtr;
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;

    DistMatrix<Real,VC,STAR> sLower(s.Grid()), zLower(z.Grid());
    SOCLowerNorms( s, sLower, orders, firstInds, cutoff );
    SOCLowerNorms( z, zLower, orders, firstInds, cutoff );

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

template<typename Real>
void ForcePairIntoSOC
(       DistMultiVec<Real>& s, 
        DistMultiVec<Real>& z,
  const DistMultiVec<Real>& w,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, 
  Real wMaxNormLimit, Int cutoff )
{
    DEBUG_ONLY(CSE cse("ForcePairIntoSOC"))

    DistMultiVec<Real> sLower(s.Comm()), zLower(z.Comm());
    SOCLowerNorms( s, sLower, orders, firstInds, cutoff );
    SOCLowerNorms( z, zLower, orders, firstInds, cutoff );

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
  template void ForcePairIntoSOC \
  (       Matrix<Real>& s, \
          Matrix<Real>& z, \
    const Matrix<Real>& w, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    Real wMaxNormLimit ); \
  template void ForcePairIntoSOC \
  (       AbstractDistMatrix<Real>& s, \
          AbstractDistMatrix<Real>& z, \
    const AbstractDistMatrix<Real>& w, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, \
    Real wMaxNormLimit, Int cutoff ); \
  template void ForcePairIntoSOC \
  (       DistMultiVec<Real>& s, \
          DistMultiVec<Real>& z, \
    const DistMultiVec<Real>& w, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    Real wMaxNormLimit, Int cutoff );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
