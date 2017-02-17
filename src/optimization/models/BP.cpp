/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include "./BP/ADMM.hpp"
#include "./BP/IPM.hpp"

namespace El {

namespace bp {

template<typename Real,
         typename=EnableIf<IsBlasScalar<Real>>>
void Helper
( const Matrix<Real>& A,
  const Matrix<Real>& b,
        Matrix<Real>& x,
  const BPCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.useIPM )
    {
        if( ctrl.useSOCP )
            bp::SOCPIPM( A, b, x, ctrl.socpIPMCtrl );
        else
            bp::LPIPM( A, b, x, ctrl.lpIPMCtrl );
    }
    else
        bp::ADMM( A, b, x, ctrl.admmCtrl );
}

template<typename Real,
         typename=DisableIf<IsBlasScalar<Real>>,
         typename=void>
void Helper
( const Matrix<Real>& A,
  const Matrix<Real>& b,
        Matrix<Real>& x,
  const BPCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.useIPM )
    {
        if( ctrl.useSOCP )
            bp::SOCPIPM( A, b, x, ctrl.socpIPMCtrl );
        else
            bp::LPIPM( A, b, x, ctrl.lpIPMCtrl );
    }
    else
        LogicError("ADMM BasisPursuit not available for this datatype");
}

template<typename Real>
void Helper
( const Matrix<Complex<Real>>& A,
  const Matrix<Complex<Real>>& b,
        Matrix<Complex<Real>>& x,
  const BPCtrl<Complex<Real>>& ctrl )
{
    EL_DEBUG_CSE
    bp::SOCPIPM( A, b, x, ctrl.ipmCtrl );
}

template<typename Real,
         typename=EnableIf<IsBlasScalar<Real>>>
void Helper
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b,
        AbstractDistMatrix<Real>& x,
  const BPCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.useIPM )
    {
        if( ctrl.useSOCP )
            bp::SOCPIPM( A, b, x, ctrl.socpIPMCtrl );
        else
            bp::LPIPM( A, b, x, ctrl.lpIPMCtrl );
    }
    else
        bp::ADMM( A, b, x, ctrl.admmCtrl );
}

template<typename Real,
         typename=DisableIf<IsBlasScalar<Real>>,
         typename=void>
void Helper
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b,
        AbstractDistMatrix<Real>& x,
  const BPCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.useIPM )
    {
        if( ctrl.useSOCP )
            bp::SOCPIPM( A, b, x, ctrl.socpIPMCtrl );
        else
            bp::LPIPM( A, b, x, ctrl.lpIPMCtrl );
    }
    else
        LogicError("ADMM BasisPursuit not available for this datatype");
}

template<typename Real>
void Helper
( const AbstractDistMatrix<Complex<Real>>& A,
  const AbstractDistMatrix<Complex<Real>>& b,
        AbstractDistMatrix<Complex<Real>>& x,
  const BPCtrl<Complex<Real>>& ctrl )
{
    EL_DEBUG_CSE
    bp::SOCPIPM( A, b, x, ctrl.ipmCtrl );
}

template<typename Real>
void Helper
( const SparseMatrix<Real>& A,
  const Matrix<Real>& b,
        Matrix<Real>& x,
  const BPCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( !ctrl.useIPM )
        LogicError("ADMM-based BP not yet supported for sparse matrices");
    if( ctrl.useSOCP )
        bp::SOCPIPM( A, b, x, ctrl.socpIPMCtrl );
    else
        bp::LPIPM( A, b, x, ctrl.lpIPMCtrl );
}

template<typename Real>
void Helper
( const SparseMatrix<Complex<Real>>& A,
  const Matrix<Complex<Real>>& b,
        Matrix<Complex<Real>>& x,
  const BPCtrl<Complex<Real>>& ctrl )
{
    EL_DEBUG_CSE
    bp::SOCPIPM( A, b, x, ctrl.ipmCtrl );
}

template<typename Real>
void Helper
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b,
        DistMultiVec<Real>& x,
  const BPCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( !ctrl.useIPM )
        LogicError("ADMM-based BP not yet supported for sparse matrices");
    if( ctrl.useSOCP )
        bp::SOCPIPM( A, b, x, ctrl.socpIPMCtrl );
    else
        bp::LPIPM( A, b, x, ctrl.lpIPMCtrl );
}

template<typename Real>
void Helper
( const DistSparseMatrix<Complex<Real>>& A,
  const DistMultiVec<Complex<Real>>& b,
        DistMultiVec<Complex<Real>>& x,
  const BPCtrl<Complex<Real>>& ctrl )
{
    EL_DEBUG_CSE
    bp::SOCPIPM( A, b, x, ctrl.ipmCtrl );
}

} // namespace bp

template<typename Field>
void BP
( const Matrix<Field>& A,
  const Matrix<Field>& b,
        Matrix<Field>& x,
  const BPCtrl<Field>& ctrl )
{ bp::Helper( A, b, x, ctrl ); }

template<typename Field>
void BP
( const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Field>& b,
        AbstractDistMatrix<Field>& x,
  const BPCtrl<Field>& ctrl )
{ bp::Helper( A, b, x, ctrl ); }

template<typename Field>
void BP
( const SparseMatrix<Field>& A,
  const Matrix<Field>& b,
        Matrix<Field>& x,
  const BPCtrl<Field>& ctrl )
{ bp::Helper( A, b, x, ctrl ); }

template<typename Field>
void BP
( const DistSparseMatrix<Field>& A,
  const DistMultiVec<Field>& b,
        DistMultiVec<Field>& x,
  const BPCtrl<Field>& ctrl )
{ bp::Helper( A, b, x, ctrl ); }

#define PROTO(Field) \
  template void BP \
  ( const Matrix<Field>& A, \
    const Matrix<Field>& b, \
          Matrix<Field>& x, \
    const BPCtrl<Field>& ctrl ); \
  template void BP \
  ( const AbstractDistMatrix<Field>& A, \
    const AbstractDistMatrix<Field>& b, \
          AbstractDistMatrix<Field>& x, \
    const BPCtrl<Field>& ctrl ); \
  template void BP \
  ( const SparseMatrix<Field>& A, \
    const Matrix<Field>& b, \
          Matrix<Field>& x, \
    const BPCtrl<Field>& ctrl ); \
  template void BP \
  ( const DistSparseMatrix<Field>& A, \
    const DistMultiVec<Field>& b, \
          DistMultiVec<Field>& x, \
    const BPCtrl<Field>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
