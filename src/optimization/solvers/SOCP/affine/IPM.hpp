/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace socp {
namespace affine {

template<typename Real>
void Mehrotra
( const Matrix<Real>& A,
  const Matrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const MehrotraCtrl<Real>& ctrl=MehrotraCtrl<Real>() );
template<typename Real>
void Mehrotra
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& b,
  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
        AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& s,
  const MehrotraCtrl<Real>& ctrl=MehrotraCtrl<Real>() );
template<typename Real>
void Mehrotra
( const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const MehrotraCtrl<Real>& ctrl=MehrotraCtrl<Real>() );
template<typename Real>
void Mehrotra
( const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
        DistMultiVec<Real>& s,
  const MehrotraCtrl<Real>& ctrl=MehrotraCtrl<Real>() );

} // namespace affine
} // namespace socp
} // namespace El
