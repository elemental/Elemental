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

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
void Helper
( const Matrix<Real>& A, 
  const Matrix<Real>& b, 
        Matrix<Real>& x,
  const BPCtrl<Real>& ctrl )
{
    DEBUG_CSE
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

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
void Helper
( const Matrix<Real>& A, 
  const Matrix<Real>& b, 
        Matrix<Real>& x,
  const BPCtrl<Real>& ctrl )
{
    DEBUG_CSE
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
    DEBUG_CSE
    bp::SOCPIPM( A, b, x, ctrl.ipmCtrl );
}

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
void Helper
( const ElementalMatrix<Real>& A, 
  const ElementalMatrix<Real>& b, 
        ElementalMatrix<Real>& x,
  const BPCtrl<Real>& ctrl )
{
    DEBUG_CSE
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

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
void Helper
( const ElementalMatrix<Real>& A, 
  const ElementalMatrix<Real>& b, 
        ElementalMatrix<Real>& x,
  const BPCtrl<Real>& ctrl )
{
    DEBUG_CSE
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
( const ElementalMatrix<Complex<Real>>& A, 
  const ElementalMatrix<Complex<Real>>& b, 
        ElementalMatrix<Complex<Real>>& x,
  const BPCtrl<Complex<Real>>& ctrl )
{
    DEBUG_CSE
    bp::SOCPIPM( A, b, x, ctrl.ipmCtrl );
}

template<typename Real>
void Helper
( const SparseMatrix<Real>& A, 
  const Matrix<Real>& b, 
        Matrix<Real>& x,
  const BPCtrl<Real>& ctrl )
{
    DEBUG_CSE
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
    DEBUG_CSE
    bp::SOCPIPM( A, b, x, ctrl.ipmCtrl );
}

template<typename Real>
void Helper
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b, 
        DistMultiVec<Real>& x,
  const BPCtrl<Real>& ctrl )
{
    DEBUG_CSE
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
    DEBUG_CSE
    bp::SOCPIPM( A, b, x, ctrl.ipmCtrl );
}

} // namespace bp

template<typename F>
void BP
( const Matrix<F>& A,
  const Matrix<F>& b,
        Matrix<F>& x,
  const BPCtrl<F>& ctrl )
{ bp::Helper( A, b, x, ctrl ); }

template<typename F>
void BP
( const ElementalMatrix<F>& A,
  const ElementalMatrix<F>& b,
        ElementalMatrix<F>& x,
  const BPCtrl<F>& ctrl )
{ bp::Helper( A, b, x, ctrl ); }

template<typename F>
void BP
( const SparseMatrix<F>& A,
  const Matrix<F>& b,
        Matrix<F>& x,
  const BPCtrl<F>& ctrl )
{ bp::Helper( A, b, x, ctrl ); }

template<typename F>
void BP
( const DistSparseMatrix<F>& A,
  const DistMultiVec<F>& b,
        DistMultiVec<F>& x,
  const BPCtrl<F>& ctrl )
{ bp::Helper( A, b, x, ctrl ); }

#define PROTO(F) \
  template void BP \
  ( const Matrix<F>& A, \
    const Matrix<F>& b, \
          Matrix<F>& x, \
    const BPCtrl<F>& ctrl ); \
  template void BP \
  ( const ElementalMatrix<F>& A, \
    const ElementalMatrix<F>& b, \
          ElementalMatrix<F>& x, \
    const BPCtrl<F>& ctrl ); \
  template void BP \
  ( const SparseMatrix<F>& A, \
    const Matrix<F>& b, \
          Matrix<F>& x, \
    const BPCtrl<F>& ctrl ); \
  template void BP \
  ( const DistSparseMatrix<F>& A, \
    const DistMultiVec<F>& b, \
          DistMultiVec<F>& x, \
    const BPCtrl<F>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
