/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_HERMITIANSOLVE_HPP
#define EL_HERMITIANSOLVE_HPP

#include EL_SYMMETRICSOLVE_INC

namespace El {

template<typename F>
inline void
HermitianSolve
( UpperOrLower uplo, Orientation orientation, Matrix<F>& A, Matrix<F>& B, 
  LDLPivotType pivotType=BUNCH_KAUFMAN_A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSolve"))
    SymmetricSolve( uplo, orientation, A, B, true, pivotType );
}

template<typename F>
inline void
HermitianSolve
( UpperOrLower uplo, Orientation orientation, 
  DistMatrix<F>& A, DistMatrix<F>& B, bool conjugate=false, 
  LDLPivotType pivotType=BUNCH_KAUFMAN_A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSolve"))
    SymmetricSolve( uplo, orientation, A, B, true, pivotType );
}

} // namespace El

#endif // ifndef EL_HERMITIANSOLVE_HPP
