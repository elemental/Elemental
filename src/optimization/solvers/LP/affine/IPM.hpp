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
void Mehrotra
( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const MehrotraCtrl<Real>& ctrl=MehrotraCtrl<Real>() );
template<typename Real>
void Mehrotra
( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
        AffineLPSolution<DistMatrix<Real>>& solution,
  const MehrotraCtrl<Real>& ctrl=MehrotraCtrl<Real>() );
template<typename Real>
void Mehrotra
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const MehrotraCtrl<Real>& ctrl=MehrotraCtrl<Real>() );
template<typename Real>
void Mehrotra
( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
        AffineLPSolution<DistMultiVec<Real>>& solution,
  const MehrotraCtrl<Real>& ctrl=MehrotraCtrl<Real>() );

} // namespace affine
} // namespace lp
} // namespace El
