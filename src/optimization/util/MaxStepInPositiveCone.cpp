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
Real MaxStepInPositiveCone
( const Matrix<Real>& s, const Matrix<Real>& ds, Real upperBound )
{
    DEBUG_ONLY(CSE cse("MaxStepInPositiveCone"))
    const Int k = s.Height();
    Real alpha = upperBound;
    for( Int i=0; i<k; ++i )
    {
        const Real si = s.Get(i,0);
        const Real dsi = ds.Get(i,0);
        if( dsi < Real(0) )
            alpha = Min(alpha,-si/dsi);
    }
    return alpha;
}

template<typename Real>
Real MaxStepInPositiveCone
( const AbstractDistMatrix<Real>& sPre, 
  const AbstractDistMatrix<Real>& dsPre, Real upperBound )
{
    DEBUG_ONLY(CSE cse("MaxStepInPositiveCone"))

    // TODO: Decide if more general intermediate distributions should be
    //       supported.
    auto sPtr = ReadProxy<Real,MC,MR>(&sPre);
    auto& s = *sPtr;

    ProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = s.ColAlign();
    control.rowAlign = s.RowAlign();
    auto dsPtr = ReadProxy<Real,MC,MR>(&dsPre,control);
    auto& ds = *dsPtr;

    Real alpha = upperBound;
    if( s.IsLocalCol(0) )
    {
        const Int kLocal = s.LocalHeight();
        for( Int iLoc=0; iLoc<kLocal; ++iLoc )
        {
            const Real si = s.GetLocal(iLoc,0);
            const Real dsi = ds.GetLocal(iLoc,0);
            if( dsi < Real(0) )
                alpha = Min(alpha,-si/dsi);
        }
    }
    return mpi::AllReduce( alpha, mpi::MIN, s.DistComm() );
}

template<typename Real>
Real MaxStepInPositiveCone
( const DistMultiVec<Real>& s, const DistMultiVec<Real>& ds, Real upperBound )
{
    DEBUG_ONLY(CSE cse("MaxStepInPositiveCone"))
    Real alpha = upperBound;
    const Int kLocal = s.LocalHeight();
    for( Int iLoc=0; iLoc<kLocal; ++iLoc )
    {
        const Real si = s.GetLocal(iLoc,0);
        const Real dsi = ds.GetLocal(iLoc,0);
        if( dsi < Real(0) )
            alpha = Min(alpha,-si/dsi);
    }
    return mpi::AllReduce( alpha, mpi::MIN, s.Comm() );
}

#define PROTO(Real) \
  template Real MaxStepInPositiveCone \
  ( const Matrix<Real>& s, const Matrix<Real>& ds, Real upperBound ); \
  template Real MaxStepInPositiveCone \
  ( const AbstractDistMatrix<Real>& s, \
    const AbstractDistMatrix<Real>& ds, Real upperBound ); \
  template Real MaxStepInPositiveCone \
  ( const DistMultiVec<Real>& s, \
    const DistMultiVec<Real>& ds, Real upperBound );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
