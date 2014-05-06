/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LDL_HPP
#define ELEM_LDL_HPP

namespace elem {
template<typename F>
void LocalLDL( DistMatrix<F,STAR,STAR>& A, bool conjugate=false );
} // namespace elem

#include "./LDL/Var3.hpp"
#include "./LDL/Pivoted.hpp"

#include "./LDL/MultiplyAfter.hpp"
#include "./LDL/SolveAfter.hpp"

#include "./LDL/Inertia.hpp"

namespace elem {

template<typename F>
inline void
LocalLDL( DistMatrix<F,STAR,STAR>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("LocalLDL"))
    ldl::Var3( A.Matrix(), conjugate );
}

template<typename F>
inline void
LDLH( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LDLH"))
    ldl::Var3( A, true );
}

template<typename F>
inline void
LDLH
( Matrix<F>& A, Matrix<F>& dSub, Matrix<Int>& pPerm, 
  LDLPivotType pivotType=BUNCH_KAUFMAN_A )
{
    DEBUG_ONLY(CallStackEntry cse("LDLH"))
    ldl::Pivoted( A, dSub, pPerm, true, pivotType );
}

template<typename F>
inline void 
LDLH( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LDLH"))
    ldl::Var3( A, true );
}

template<typename F,Dist UPerm>
inline void
LDLH
( DistMatrix<F>& A, 
  DistMatrix<F,MD,STAR>& dSub, 
  DistMatrix<Int,UPerm,STAR>& pPerm,
  LDLPivotType pivotType=BUNCH_KAUFMAN_A )
{
    DEBUG_ONLY(CallStackEntry cse("LDLH"))
    ldl::Pivoted( A, dSub, pPerm, true, pivotType );
}

template<typename F>
inline void
LDLT( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LDLT"))
    ldl::Var3( A, false );
}

template<typename F>
inline void
LDLT
( Matrix<F>& A, Matrix<F>& dSub, Matrix<Int>& pPerm, 
  LDLPivotType pivotType=BUNCH_KAUFMAN_A )
{
    DEBUG_ONLY(CallStackEntry cse("LDLT"))
    ldl::Pivoted( A, dSub, pPerm, false, pivotType );
}

template<typename F>
inline void 
LDLT( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LDLT"))
    ldl::Var3( A, false );
}

template<typename F,Dist UPerm>
inline void
LDLT
( DistMatrix<F>& A, 
  DistMatrix<F,MD,STAR>& dSub, 
  DistMatrix<Int,UPerm,STAR>& pPerm,
  LDLPivotType pivotType=BUNCH_KAUFMAN_A )
{
    DEBUG_ONLY(CallStackEntry cse("LDLT"))
    ldl::Pivoted( A, dSub, pPerm, false, pivotType );
}

} // namespace elem

#endif // ifndef ELEM_LDL_HPP
