/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./BidiagSVD/QR.hpp"

namespace El {

namespace bidiag_svd {

template<typename Real>
Real APosterioriThreshold
( Int m, Int n,
  const Real& twoNorm,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE
    Real thresh = Max(m,n)*twoNorm*limits::Epsilon<Real>();
    auto tol = ctrl.tol;
    auto tolType = ctrl.tolType;
    if( tolType == ABSOLUTE_SING_VAL_TOL )
    {
        if( tol == Real(0) )
        {
            // Nothing was specified...keep the default relative tolerance
        }
        else
        {
            thresh = tol;
        }
    }
    else if( tolType == RELATIVE_TO_MAX_SING_VAL_TOL )
    {
        if( tol == Real(0) )
        {
            // Nothing was specified...keep the default relative tolerance
        }
        else
        {
            thresh = tol*twoNorm;
        }
    }
    else
    {
        // Keeping all singular values relatively accurate means that none
        // can be truncated
        thresh = 0;
    }
    return thresh;
}

} // namespace bidiag_svd

// TODO(poulson): Generalize to complex bidiagonal via unitary diagonal
// rotations
template<typename Real,typename>
BidiagSVDInfo
BidiagSVD
( Matrix<Real>& mainDiag,
  Matrix<Real>& superDiag,
  Matrix<Real>& s,
  const BidiagSVDCtrl<Real>& ctrl )
{
    DEBUG_CSE
    BidiagSVDInfo info;
    s = mainDiag;
    info.qrInfo = bidiag_svd::QRAlg( s, superDiag, ctrl );
    Sort( mainDiag, ASCENDING );
    return info;
}

template<typename F>
BidiagSVDInfo
BidiagSVD
( Matrix<Base<F>>& mainDiag,
  Matrix<Base<F>>& superDiag,
  Matrix<F>& U,
  Matrix<Base<F>>& s,
  Matrix<F>& V,
  const BidiagSVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    const Int n = mainDiag.Height();
    BidiagSVDInfo info;

    s = mainDiag;
    info.qrInfo = bidiag_svd::QRAlg( s, superDiag, U, V, ctrl );
    if( ctrl.wantU || ctrl.wantV )
    {
        auto sortPairs = TaggedSort( mainDiag, ASCENDING );
        for( Int j=0; j<n; ++j )
            mainDiag(j) = sortPairs[j].value;
        if( ctrl.wantU )
            ApplyTaggedSortToEachRow( sortPairs, U );
        if( ctrl.wantV )
            ApplyTaggedSortToEachRow( sortPairs, V );
    }
    else
    {
        Sort( s, ASCENDING );
    }

    return info;
}

#define PROTO(F) \
  template BidiagSVDInfo BidiagSVD \
  ( Matrix<Base<F>>& mainDiag, \
    Matrix<Base<F>>& superDiag, \
    Matrix<F>& U, \
    Matrix<Base<F>>& s, \
    Matrix<F>& V, \
    const BidiagSVDCtrl<Base<F>>& ctrl ); \
  template bidiag_svd::QRInfo bidiag_svd::QRAlg \
  ( Matrix<Base<F>>& mainDiag, \
    Matrix<Base<F>>& superDiag, \
    Matrix<F>& U, \
    Matrix<F>& V, \
    const BidiagSVDCtrl<Base<F>>& ctrl );

#define PROTO_REAL(Real) \
  PROTO(Real) \
  template BidiagSVDInfo BidiagSVD \
  ( Matrix<Real>& mainDiag, \
    Matrix<Real>& superDiag, \
    Matrix<Real>& s, \
    const BidiagSVDCtrl<Real>& ctrl ); \
  template bidiag_svd::QRInfo bidiag_svd::QRAlg \
  ( Matrix<Real>& mainDiag, \
    Matrix<Real>& superDiag, \
    const BidiagSVDCtrl<Real>& ctrl ); \
  template Real bidiag_svd::APosterioriThreshold \
  ( Int m, Int n, \
    const Real& twoNorm, \
    const BidiagSVDCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
