/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_HPDSOLVE_HPP
#define EL_HPDSOLVE_HPP

#include EL_TRSM_INC
#include EL_CHOLESKY_INC

namespace El {

template<typename F>
inline void
HPDSolve
( UpperOrLower uplo, Orientation orientation, Matrix<F>& A, Matrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("HPDSolve"))
    Cholesky( uplo, A );
    cholesky::SolveAfter( uplo, orientation, A, B );
}

template<typename F>
inline void
HPDSolve
( UpperOrLower uplo, Orientation orientation, 
  DistMatrix<F>& A, DistMatrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("HPDSolve"))
    Cholesky( uplo, A );
    cholesky::SolveAfter( uplo, orientation, A, B );
}

} // namespace El

#endif // ifndef EL_HPDSOLVE_HPP
