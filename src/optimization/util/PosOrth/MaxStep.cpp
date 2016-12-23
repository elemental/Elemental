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
Real MaxStep
( const Matrix<Real>& s,
  const Matrix<Real>& ds,
        Real upperBound )
{
    EL_DEBUG_CSE
    const Int k = s.Height();
    Real alpha = upperBound;
    for( Int i=0; i<k; ++i )
    {
        if( ds(i) < Real(0) )
            alpha = Min(alpha,-s(i)/ds(i));
    }
    return alpha;
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Real MaxStep
( const AbstractDistMatrix<Real>& sPre,
  const AbstractDistMatrix<Real>& dsPre,
  Real upperBound )
{
    EL_DEBUG_CSE

    // TODO(poulson): Decide if more general intermediate distributions should
    // be supported.
    DistMatrixReadProxy<Real,Real,MC,MR> sProx( sPre );
    auto& s = sProx.GetLocked();

    ElementalProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = s.ColAlign();
    control.rowAlign = s.RowAlign();

    DistMatrixReadProxy<Real,Real,MC,MR> dsProx( dsPre, control );
    auto& ds = dsProx.GetLocked();

    const Real* sBuf = s.LockedBuffer();
    const Real* dsBuf = ds.LockedBuffer();

    Real alpha = upperBound;
    if( s.IsLocalCol(0) )
    {
        const Int kLocal = s.LocalHeight();
        for( Int iLoc=0; iLoc<kLocal; ++iLoc )
        {
            const Real si = sBuf[iLoc];
            const Real dsi = dsBuf[iLoc];
            if( dsi < Real(0) )
                alpha = Min(alpha,-si/dsi);
        }
    }
    return mpi::AllReduce( alpha, mpi::MIN, s.DistComm() );
}

template<typename Real,
         typename/*=EnableIf<IsReal<Real>>*/>
Real MaxStep
( const DistMultiVec<Real>& s,
  const DistMultiVec<Real>& ds,
  Real upperBound )
{
    EL_DEBUG_CSE
    const Int kLocal = s.LocalHeight();
    const Real* sBuf = s.LockedMatrix().LockedBuffer();
    const Real* dsBuf = ds.LockedMatrix().LockedBuffer();

    Real alpha = upperBound;
    for( Int iLoc=0; iLoc<kLocal; ++iLoc )
    {
        const Real si = sBuf[iLoc];
        const Real dsi = dsBuf[iLoc];
        if( dsi < Real(0) )
            alpha = Min(alpha,-si/dsi);
    }
    return mpi::AllReduce( alpha, mpi::MIN, s.Grid().Comm() );
}

#define PROTO(Real) \
  template Real MaxStep \
  ( const Matrix<Real>& s, \
    const Matrix<Real>& ds, \
    Real upperBound ); \
  template Real MaxStep \
  ( const AbstractDistMatrix<Real>& s, \
    const AbstractDistMatrix<Real>& ds, \
    Real upperBound ); \
  template Real MaxStep \
  ( const DistMultiVec<Real>& s, \
    const DistMultiVec<Real>& ds, \
    Real upperBound );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace pos_orth
} // namespace El
