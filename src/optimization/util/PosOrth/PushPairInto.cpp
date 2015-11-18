/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace pos_orth {

template<typename Real,typename>
void PushPairInto
(       Matrix<Real>& s, 
        Matrix<Real>& z,
  const Matrix<Real>& w,
  Real wMaxNormLimit )
{
    DEBUG_ONLY(CSE cse("pos_orth::PushPairInto"))
    const Int height = s.Height();
    const Real maxMod = Pow(Epsilon<Real>(),Real(0.5));
    for( Int i=0; i<height; ++i )
    {
        if( w.Get(i,0) > wMaxNormLimit )
        {
            // TODO: Switch to a non-adhoc modification     
            z.Update( i, 0, Min(Real(1)/wMaxNormLimit,maxMod) );
        }
    }
}

template<typename Real,typename>
void PushPairInto
(       ElementalMatrix<Real>& sPre, 
        ElementalMatrix<Real>& zPre,
  const ElementalMatrix<Real>& wPre,
  Real wMaxNormLimit )
{
    DEBUG_ONLY(CSE cse("pos_orth::PushPairInto"))
    AssertSameGrids( sPre, zPre, wPre );
    const Real maxMod = Pow(Epsilon<Real>(),Real(0.5));

    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;

    auto sPtr = ReadWriteProxy<Real,VC,STAR>(&sPre,ctrl); 
    auto zPtr = ReadWriteProxy<Real,VC,STAR>(&zPre,ctrl); 
    auto wPtr = ReadProxy<Real,VC,STAR>(&wPre,ctrl);
    auto& s = *sPtr;
    auto& z = *zPtr;
    auto& w = *wPtr;

    const Int localHeight = s.LocalHeight();
    const Real* wBuf = w.LockedBuffer();
    Real* zBuf = z.Buffer();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        if( wBuf[iLoc] > wMaxNormLimit )
        {
            // TODO: Switch to a non-adhoc modification     
            zBuf[iLoc] += Min(Real(1)/wMaxNormLimit,maxMod);
        }
    }
}

template<typename Real,typename>
void PushPairInto
(       DistMultiVec<Real>& s, 
        DistMultiVec<Real>& z,
  const DistMultiVec<Real>& w,
  Real wMaxNormLimit )
{
    DEBUG_ONLY(CSE cse("pos_orth::PushPairInto"))
    const Real maxMod = Pow(Epsilon<Real>(),Real(0.5));
    const int localHeight = s.LocalHeight();
    const Real* wBuf = w.LockedMatrix().LockedBuffer();
    Real* zBuf = z.Matrix().Buffer();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        if( wBuf[iLoc] > wMaxNormLimit )
        {
            // TODO: Switch to a non-adhoc modification     
            zBuf[iLoc] += Min(Real(1)/wMaxNormLimit,maxMod);
        }
    }
}

#define PROTO(Real) \
  template void PushPairInto \
  (       Matrix<Real>& s, \
          Matrix<Real>& z, \
    const Matrix<Real>& w, \
    Real wMaxNormLimit ); \
  template void PushPairInto \
  (       ElementalMatrix<Real>& s, \
          ElementalMatrix<Real>& z, \
    const ElementalMatrix<Real>& w, \
    Real wMaxNormLimit ); \
  template void PushPairInto \
  (       DistMultiVec<Real>& s, \
          DistMultiVec<Real>& z, \
    const DistMultiVec<Real>& w, \
    Real wMaxNormLimit );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace pos_orth
} // namespace El
