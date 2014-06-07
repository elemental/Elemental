/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include "./LDL/Var3.hpp"
#include "./LDL/Pivoted.hpp"

#include "./LDL/MultiplyAfter.hpp"
#include "./LDL/SolveAfter.hpp"

#include "./LDL/Inertia.hpp"

namespace El {

template<typename F>
void LDL( Matrix<F>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("LDL"))
    ldl::Var3( A, conjugate );
}

template<typename F>
void LDL( DistMatrix<F>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("LDL"))
    ldl::Var3( A, conjugate );
}

template<typename F>
void LDL
( Matrix<F>& A, Matrix<F>& dSub, 
  Matrix<Int>& pPerm, bool conjugate, LDLPivotType pivotType )
{
    DEBUG_ONLY(CallStackEntry cse("LDL"))
    ldl::Pivoted( A, dSub, pPerm, conjugate, pivotType );
}

template<typename F,Dist UPerm>
void LDL
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& dSub, 
  DistMatrix<Int,UPerm,STAR>& pPerm, bool conjugate, LDLPivotType pivotType )
{
    DEBUG_ONLY(CallStackEntry cse("LDL"))
    ldl::Pivoted( A, dSub, pPerm, conjugate, pivotType );
}

#define PROTO(F) \
  template void LDL( Matrix<F>& A, bool conjugate ); \
  template void LDL( DistMatrix<F>& A, bool conjugate ); \
  template void LDL \
  ( Matrix<F>& A, Matrix<F>& dSub, \
    Matrix<Int>& pPerm, bool conjugate, LDLPivotType pivotType ); \
  template void LDL \
  ( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& dSub, \
    DistMatrix<Int,VC,STAR>& pPerm, bool conjugate, LDLPivotType pivotType ); \
  template InertiaType ldl::Inertia \
  ( const Matrix<Base<F>>& d, const Matrix<F>& dSub ); \
  template InertiaType ldl::Inertia \
  ( const DistMatrix<Base<F>,MD,STAR>& d, const DistMatrix<F,MD,STAR>& dSub ); \
  template void ldl::MultiplyAfter \
  ( const Matrix<F>& A, Matrix<F>& B, bool conjugated ); \
  template void ldl::MultiplyAfter \
  ( const DistMatrix<F>& A, DistMatrix<F>& B, bool conjugated ); \
  template void ldl::MultiplyAfter \
  ( const Matrix<F>& A, const Matrix<F>& dSub, \
    const Matrix<Int>& pPerm, Matrix<F>& B, bool conjugated ); \
  template void ldl::MultiplyAfter \
  ( const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& dSub, \
    const DistMatrix<Int,VC,STAR>& pPerm, DistMatrix<F>& B, bool conjugated ); \
  template void ldl::SolveAfter \
  ( const Matrix<F>& A, Matrix<F>& B, bool conjugated ); \
  template void ldl::SolveAfter \
  ( const DistMatrix<F>& A, DistMatrix<F>& B, bool conjugated ); \
  template void ldl::SolveAfter \
  ( const Matrix<F>& A, const Matrix<F>& dSub, \
    const Matrix<Int>& pPerm, Matrix<F>& B, bool conjugated ); \
  template void ldl::SolveAfter \
  ( const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& dSub, \
    const DistMatrix<Int,VC,STAR>& pPerm, DistMatrix<F>& B, bool conjugated );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
