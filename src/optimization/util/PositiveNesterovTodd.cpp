/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Find the Nesterov-Todd scaling point w such that 
//
//   Q_w z = diag(w)^2 z = s,
//
// which implies that w_i = sqrt(s_i/z_i).

template<typename Real>
void PositiveNesterovTodd
( const Matrix<Real>& s, 
  const Matrix<Real>& z,
        Matrix<Real>& w )
{
    DEBUG_ONLY(CSE cse("PositiveNesterovTodd"))
    const Int k = s.Height();
    w.Resize( k, 1 );
    for( Int i=0; i<k; ++i )
        w.Set( i, 0, Sqrt(s.Get(i,0)/z.Get(i,0)) );
}

template<typename Real>
void PositiveNesterovTodd
( const AbstractDistMatrix<Real>& sPre, 
  const AbstractDistMatrix<Real>& zPre,
        AbstractDistMatrix<Real>& wPre )
{
    DEBUG_ONLY(CSE cse("PositiveNesterovTodd"))
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
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        w.SetLocal( iLoc, 0, Sqrt(s.GetLocal(iLoc,0)/z.GetLocal(iLoc,0)) );
}

template<typename Real>
void PositiveNesterovTodd
( const DistMultiVec<Real>& s, 
  const DistMultiVec<Real>& z,
        DistMultiVec<Real>& w )
{
    DEBUG_ONLY(CSE cse("PositiveNesterovTodd"))
    w.SetComm( s.Comm() );
    w.Resize( s.Height(), 1 );
    const Int localHeight = w.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        w.SetLocal( iLoc, 0, Sqrt(s.GetLocal(iLoc,0)/z.GetLocal(iLoc,0)) );
}

#define PROTO(Real) \
  template void PositiveNesterovTodd \
  ( const Matrix<Real>& s, \
    const Matrix<Real>& z, \
          Matrix<Real>& w ); \
  template void PositiveNesterovTodd \
  ( const AbstractDistMatrix<Real>& s, \
    const AbstractDistMatrix<Real>& z, \
          AbstractDistMatrix<Real>& w ); \
  template void PositiveNesterovTodd \
  ( const DistMultiVec<Real>& s, \
    const DistMultiVec<Real>& z, \
          DistMultiVec<Real>& w );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
