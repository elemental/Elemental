/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace lp {
namespace affine {

template<typename Real>
IPMInfo<Real> IPM
( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl=IPMCtrl<Real>() );
template<typename Real>
IPMInfo<Real> IPM
( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
        AffineLPSolution<DistMatrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl=IPMCtrl<Real>() );
template<typename Real>
IPMInfo<Real> IPM
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl=IPMCtrl<Real>() );
template<typename Real>
IPMInfo<Real> IPM
( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
        AffineLPSolution<DistMultiVec<Real>>& solution,
  const IPMCtrl<Real>& ctrl=IPMCtrl<Real>() );

} // namespace affine
} // namespace lp
} // namespace El
