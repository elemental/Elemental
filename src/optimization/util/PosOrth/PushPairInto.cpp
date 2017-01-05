/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace pos_orth {

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void PushPairInto
(       Matrix<Real>& s,
        Matrix<Real>& z,
  const Matrix<Real>& w,
  Real wRatioLimit )
{
    EL_DEBUG_CSE
    const Int height = s.Height();
    const Real ratioLimit = wRatioLimit;
    for( Int i=0; i<height; ++i )
    {
        if( s(i) / z(i) > ratioLimit )
        {
            z(i) = s(i) / ratioLimit;
        }
        else if( z(i) / s(i) > ratioLimit )
        {
            s(i) = z(i) / ratioLimit;
        }
    }
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void PushPairInto
(       AbstractDistMatrix<Real>& sPre,
        AbstractDistMatrix<Real>& zPre,
  const AbstractDistMatrix<Real>& wPre,
  Real wRatioLimit )
{
    EL_DEBUG_CSE
    AssertSameGrids( sPre, zPre, wPre );
    const Real ratioLimit = wRatioLimit;

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    DistMatrixWriteProxy<Real,Real,VC,STAR>
      sProx( sPre, ctrl ),
      zProx( zPre, ctrl );
    /*
    DistMatrixReadProxy<Real,Real,VC,STAR>
      wProx( wPre, ctrl );
    */
    auto& s = sProx.Get();
    auto& z = zProx.Get();
    //auto& w = wProx.GetLocked();

    const Int localHeight = s.LocalHeight();
    //const Real* wBuf = w.LockedBuffer();
    Real* sBuf = s.Buffer();
    Real* zBuf = z.Buffer();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        if( sBuf[iLoc] / zBuf[iLoc] > ratioLimit )
        {
            zBuf[iLoc] = sBuf[iLoc] / ratioLimit; 
        }
        else if( zBuf[iLoc] / sBuf[iLoc] > ratioLimit )
        {
            sBuf[iLoc] = zBuf[iLoc] / ratioLimit;
        }
    }
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
void PushPairInto
(       DistMultiVec<Real>& s,
        DistMultiVec<Real>& z,
  const DistMultiVec<Real>& w,
  Real wRatioLimit )
{
    EL_DEBUG_CSE
    const Real ratioLimit = wRatioLimit;

    const int localHeight = s.LocalHeight();
    //const Real* wBuf = w.LockedMatrix().LockedBuffer();
    Real* sBuf = s.Matrix().Buffer();
    Real* zBuf = z.Matrix().Buffer();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        if( sBuf[iLoc] / zBuf[iLoc] > ratioLimit )
        {
            zBuf[iLoc] = sBuf[iLoc] / ratioLimit; 
        }
        else if( zBuf[iLoc] / sBuf[iLoc] > ratioLimit )
        {
            sBuf[iLoc] = zBuf[iLoc] / ratioLimit;
        }
    }
}

#define PROTO(Real) \
  template void PushPairInto \
  (       Matrix<Real>& s, \
          Matrix<Real>& z, \
    const Matrix<Real>& w, \
    Real wRatioLimit ); \
  template void PushPairInto \
  (       AbstractDistMatrix<Real>& s, \
          AbstractDistMatrix<Real>& z, \
    const AbstractDistMatrix<Real>& w, \
    Real wRatioLimit ); \
  template void PushPairInto \
  (       DistMultiVec<Real>& s, \
          DistMultiVec<Real>& z, \
    const DistMultiVec<Real>& w, \
    Real wRatioLimit );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace pos_orth
} // namespace El
