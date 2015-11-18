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

// Find the Nesterov-Todd scaling point w such that 
//
//   Q_w z = diag(w)^2 z = s,
//
// which implies that w_i = sqrt(s_i/z_i).

template<typename Real,typename>
void NesterovTodd
( const Matrix<Real>& s, 
  const Matrix<Real>& z,
        Matrix<Real>& w )
{
    DEBUG_ONLY(CSE cse("pos_orth::NesterovTodd"))
    const Int k = s.Height();
    w.Resize( k, 1 );
    const Real* sBuf = s.LockedBuffer();
    const Real* zBuf = z.LockedBuffer();
          Real* wBuf = w.Buffer();

    for( Int i=0; i<k; ++i )
        wBuf[i] = Sqrt(sBuf[i]/zBuf[i]);
}

template<typename Real,typename>
void NesterovTodd
( const ElementalMatrix<Real>& sPre, 
  const ElementalMatrix<Real>& zPre,
        ElementalMatrix<Real>& wPre )
{
    DEBUG_ONLY(CSE cse("pos_orth::NesterovTodd"))
    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;
    auto sPtr = ReadProxy<Real,VC,STAR>(&sPre,ctrl); 
    auto zPtr = ReadProxy<Real,VC,STAR>(&zPre,ctrl);
    auto wPtr = WriteProxy<Real,VC,STAR>(&wPre,ctrl);
    auto& s = *sPtr;
    auto& z = *zPtr;
    auto& w = *wPtr;

    w.Resize( s.Height(), 1 );
    const Int localHeight = w.LocalHeight();
    const Real* sBuf = s.LockedBuffer();
    const Real* zBuf = z.LockedBuffer();
          Real* wBuf = w.Buffer();

    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        wBuf[iLoc] = Sqrt(sBuf[iLoc]/zBuf[iLoc]);
}

template<typename Real,typename>
void NesterovTodd
( const DistMultiVec<Real>& s, 
  const DistMultiVec<Real>& z,
        DistMultiVec<Real>& w )
{
    DEBUG_ONLY(CSE cse("pos_orth::NesterovTodd"))
    w.SetComm( s.Comm() );
    w.Resize( s.Height(), 1 );
    const Real* sBuf = s.LockedMatrix().LockedBuffer();
    const Real* zBuf = z.LockedMatrix().LockedBuffer();
          Real* wBuf = w.Matrix().Buffer();

    const Int localHeight = w.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        wBuf[iLoc] = Sqrt(sBuf[iLoc]/zBuf[iLoc]);
}

#define PROTO(Real) \
  template void NesterovTodd \
  ( const Matrix<Real>& s, \
    const Matrix<Real>& z, \
          Matrix<Real>& w ); \
  template void NesterovTodd \
  ( const ElementalMatrix<Real>& s, \
    const ElementalMatrix<Real>& z, \
          ElementalMatrix<Real>& w ); \
  template void NesterovTodd \
  ( const DistMultiVec<Real>& s, \
    const DistMultiVec<Real>& z, \
          DistMultiVec<Real>& w );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace pos_orth
} // namespace El
