/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./LDL/Var3.hpp"

#include "./LDL/Pivoted.hpp"

#include "./LDL/MultiplyAfter.hpp"
#include "./LDL/SolveAfter.hpp"

#include "./LDL/Inertia.hpp"

namespace El {

// Unpivoted
// =========

template<typename F>
void LDL( Matrix<F>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("LDL"))
    ldl::Var3( A, conjugate );
}

template<typename F>
void LDL( AbstractDistMatrix<F>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("LDL"))
    ldl::Var3( A, conjugate );
}

// Pivoted
// =======

template<typename F>
void LDL
( Matrix<F>& A, Matrix<F>& dSub, 
  Matrix<Int>& p, bool conjugate, const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LDL"))
    ldl::Pivoted( A, dSub, p, conjugate, ctrl );
}

template<typename F>
void LDL
( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& dSub, 
  AbstractDistMatrix<Int>& p, bool conjugate, 
  const LDLPivotCtrl<Base<F>>& ctrl ) 
{
    DEBUG_ONLY(CallStackEntry cse("LDL"))
    ldl::Pivoted( A, dSub, p, conjugate, ctrl );
}

#define PROTO(F) \
  template void LDL( Matrix<F>& A, bool conjugate ); \
  template void LDL( AbstractDistMatrix<F>& A, bool conjugate ); \
  template void LDL \
  ( Matrix<F>& A, Matrix<F>& dSub, \
    Matrix<Int>& p, bool conjugate, const LDLPivotCtrl<Base<F>>& ctrl ); \
  template void LDL \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& dSub, \
    AbstractDistMatrix<Int>& p, bool conjugate, \
    const LDLPivotCtrl<Base<F>>& ctrl ); \
  template InertiaType ldl::Inertia \
  ( const Matrix<Base<F>>& d, const Matrix<F>& dSub ); \
  template InertiaType ldl::Inertia \
  ( const AbstractDistMatrix<Base<F>>& d, const AbstractDistMatrix<F>& dSub ); \
  template void ldl::MultiplyAfter \
  ( const Matrix<F>& A, Matrix<F>& B, bool conjugated ); \
  template void ldl::MultiplyAfter \
  ( const AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B, \
    bool conjugated ); \
  template void ldl::MultiplyAfter \
  ( const Matrix<F>& A, const Matrix<F>& dSub, \
    const Matrix<Int>& p, Matrix<F>& B, bool conjugated ); \
  template void ldl::MultiplyAfter \
  ( const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& dSub, \
    const AbstractDistMatrix<Int>& p, AbstractDistMatrix<F>& B, \
    bool conjugated ); \
  template void ldl::SolveAfter \
  ( const Matrix<F>& A, Matrix<F>& B, bool conjugated ); \
  template void ldl::SolveAfter \
  ( const AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B, \
    bool conjugated ); \
  template void ldl::SolveAfter \
  ( const Matrix<F>& A, const Matrix<F>& dSub, \
    const Matrix<Int>& p, Matrix<F>& B, bool conjugated ); \
  template void ldl::SolveAfter \
  ( const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& dSub, \
    const AbstractDistMatrix<Int>& p, AbstractDistMatrix<F>& B, \
     bool conjugated );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
