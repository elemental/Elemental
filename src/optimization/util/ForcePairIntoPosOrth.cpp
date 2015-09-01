/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename Real>
void ForcePairIntoPosOrth
(       Matrix<Real>& s, 
        Matrix<Real>& z,
  const Matrix<Real>& w,
  Real wMaxNormLimit )
{
    DEBUG_ONLY(CSE cse("ForcePairIntoPosOrth"))
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

template<typename Real>
void ForcePairIntoPosOrth
(       ElementalMatrix<Real>& sPre, 
        ElementalMatrix<Real>& zPre,
  const ElementalMatrix<Real>& wPre,
  Real wMaxNormLimit )
{
    DEBUG_ONLY(CSE cse("ForcePairIntoPosOrth"))
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
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        if( w.GetLocal(iLoc,0) > wMaxNormLimit )
        {
            // TODO: Switch to a non-adhoc modification     
            z.UpdateLocal( iLoc, 0, Min(Real(1)/wMaxNormLimit,maxMod) );
        }
    }
}

template<typename Real>
void ForcePairIntoPosOrth
(       DistMultiVec<Real>& s, 
        DistMultiVec<Real>& z,
  const DistMultiVec<Real>& w,
  Real wMaxNormLimit )
{
    DEBUG_ONLY(CSE cse("ForcePairIntoPosOrth"))
    const Real maxMod = Pow(Epsilon<Real>(),Real(0.5));
    const int localHeight = s.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        if( w.GetLocal(iLoc,0) > wMaxNormLimit )
        {
            // TODO: Switch to a non-adhoc modification     
            z.UpdateLocal( iLoc, 0, Min(Real(1)/wMaxNormLimit,maxMod) );
        }
    }
}

#define PROTO(Real) \
  template void ForcePairIntoPosOrth \
  (       Matrix<Real>& s, \
          Matrix<Real>& z, \
    const Matrix<Real>& w, \
    Real wMaxNormLimit ); \
  template void ForcePairIntoPosOrth \
  (       ElementalMatrix<Real>& s, \
          ElementalMatrix<Real>& z, \
    const ElementalMatrix<Real>& w, \
    Real wMaxNormLimit ); \
  template void ForcePairIntoPosOrth \
  (       DistMultiVec<Real>& s, \
          DistMultiVec<Real>& z, \
    const DistMultiVec<Real>& w, \
    Real wMaxNormLimit );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
