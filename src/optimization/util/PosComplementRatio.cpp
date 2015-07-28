/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Compute Max( s o z ) / Min( s o z ) to determine if we need to recenter

template<typename Real>
Real PosComplementRatio
( const Matrix<Real>& s, const Matrix<Real>& z )
{
    DEBUG_ONLY(CSE cse("PosComplementRatio"))
    const Int k = s.Height();

    Real maxProd = 0;
    for( Int i=0; i<k; ++i )
        maxProd = Max( s.Get(i,0)*z.Get(i,0), maxProd );

    Real minProd = maxProd;
    for( Int i=0; i<k; ++i )
        minProd = Min( s.Get(i,0)*z.Get(i,0), minProd );

    return maxProd/minProd;
}

template<typename Real>
Real PosComplementRatio
( const AbstractDistMatrix<Real>& sPre, const AbstractDistMatrix<Real>& zPre )
{
    DEBUG_ONLY(CSE cse("PosComplementRatio"))
    ProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = 0;
    auto sPtr = ReadProxy<Real,VC,STAR>(&sPre,ctrl); 
    auto zPtr = ReadProxy<Real,VC,STAR>(&zPre,ctrl);
    auto& s = *sPtr;
    auto& z = *zPtr;

    const Int localHeight = s.LocalHeight();

    Real maxLocProd = 0;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        maxLocProd = Max( s.GetLocal(iLoc,0)*z.GetLocal(iLoc,0), maxLocProd );
    const Real maxProd = mpi::AllReduce( maxLocProd, mpi::MAX, s.DistComm() );
    
    Real minLocProd = maxProd;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        minLocProd = Min( s.GetLocal(iLoc,0)*z.GetLocal(iLoc,0), minLocProd );
    const Real minProd = mpi::AllReduce( minLocProd, mpi::MIN, s.DistComm() );

    return maxProd/minProd;
}

template<typename Real>
Real PosComplementRatio
( const DistMultiVec<Real>& s, const DistMultiVec<Real>& z )
{
    DEBUG_ONLY(CSE cse("PosComplementRatio"))
    const Int localHeight = s.LocalHeight();

    Real maxLocProd = 0;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        maxLocProd = Max( s.GetLocal(iLoc,0)*z.GetLocal(iLoc,0), maxLocProd );
    const Real maxProd = mpi::AllReduce( maxLocProd, mpi::MAX, s.Comm() );
    
    Real minLocProd = maxProd;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        minLocProd = Min( s.GetLocal(iLoc,0)*z.GetLocal(iLoc,0), minLocProd );
    const Real minProd = mpi::AllReduce( minLocProd, mpi::MIN, s.Comm() );

    return maxProd/minProd;
}

#define PROTO(Real) \
  template Real PosComplementRatio \
  ( const Matrix<Real>& s, const Matrix<Real>& z ); \
  template Real PosComplementRatio \
  ( const AbstractDistMatrix<Real>& s, const AbstractDistMatrix<Real>& z ); \
  template Real PosComplementRatio \
  ( const DistMultiVec<Real>& s, const DistMultiVec<Real>& z );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
