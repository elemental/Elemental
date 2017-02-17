/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_PSEUDOSPECTRA_ANALYTIC_HPP
#define EL_PSEUDOSPECTRA_ANALYTIC_HPP

#include "./Util.hpp"

namespace El {
namespace pspec {

template<typename Real>
void Analytic
( const Matrix<Complex<Real>>& w,
  const Matrix<Complex<Real>>& shifts,
        Matrix<Real         >& invNorms,
        SnapshotCtrl& snapCtrl )
{
    EL_DEBUG_CSE
    using namespace pspec;
    typedef Complex<Real> C;
    const Int n = w.Height();
    const Int numShifts = shifts.Height();
    const Real normCap = NormCap<Real>();

    Zeros( invNorms, numShifts, 1 );
    if( n == 0 )
        return;

    for( Int j=0; j<numShifts; ++j )
    {
        const C shift = shifts(j);
        Real minDist = Abs(shift-w(0));
        for( Int k=1; k<n; ++k )
        {
            const Real dist = Abs(shift-w(k));
            minDist = Min(dist,minDist);
        }
        Real alpha = Real(1)/minDist;
        if( !limits::IsFinite(alpha) || alpha >= normCap )
            alpha = normCap;
        invNorms(j) = alpha;
    }

    snapCtrl.itCounts = false;
    Matrix<Int> itCounts;
    FinalSnapshot( invNorms, itCounts, snapCtrl );
}

template<typename Real>
void Analytic
( const AbstractDistMatrix<Complex<Real>>& w,
  const AbstractDistMatrix<Complex<Real>>& shiftsPre,
        AbstractDistMatrix<Real>& invNormsPre,
        SnapshotCtrl& snapCtrl )
{
    EL_DEBUG_CSE
    using namespace pspec;
    typedef Complex<Real> C;

    DistMatrixReadProxy<C,C,VR,STAR> shiftsProx( shiftsPre );
    auto& shifts = shiftsProx.GetLocked();
    auto& shiftsLoc = shifts.LockedMatrix();

    ElementalProxyCtrl ctrl;
    ctrl.colConstrain = true;
    ctrl.colAlign = shifts.ColAlign();

    DistMatrixWriteProxy<Real,Real,VR,STAR> invNormsProx( invNormsPre, ctrl );
    auto& invNorms = invNormsProx.Get();

    const Int n = w.Height();
    const Int numShifts = shifts.Height();
    const Real normCap = NormCap<Real>();
    const Grid& g = w.Grid();

    invNorms.AlignWith( shifts );
    Zeros( invNorms, numShifts, 1 );
    if( n == 0 )
        return;
    auto& invNormsLoc = invNorms.Matrix();

    DistMatrix<C,STAR,STAR> w_STAR_STAR( w );
    auto& wLoc = w_STAR_STAR.LockedMatrix();

    const Int numLocShifts = shifts.LocalHeight();
    for( Int jLoc=0; jLoc<numLocShifts; ++jLoc )
    {
        const C shift = shiftsLoc(jLoc);
        Real minDist = Abs(shift-wLoc(0));
        for( Int k=1; k<n; ++k )
        {
            const Real dist = Abs(shift-wLoc(k));
            minDist = Min(dist,minDist);
        }
        Real alpha = Real(1)/minDist;
        if( !limits::IsFinite(alpha) || alpha >= normCap )
            alpha = normCap;
        invNormsLoc(jLoc) = alpha;
    }

    snapCtrl.itCounts = false;
    DistMatrix<Int,VR,STAR> itCounts(g);
    FinalSnapshot( invNorms, itCounts, snapCtrl );
}

} // namespace pspec
} // namespace El

#endif // ifndef EL_PSEUDOSPECTRA_ANALYTIC_HPP
