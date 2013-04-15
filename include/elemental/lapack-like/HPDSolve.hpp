/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HPDSOLVE_HPP
#define LAPACK_HPDSOLVE_HPP

#include "elemental/blas-like/level3/Trsm.hpp"
#include "elemental/lapack-like/Cholesky.hpp"

namespace elem {

template<typename F>
inline void
HPDSolve
( UpperOrLower uplo, Orientation orientation, Matrix<F>& A, Matrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("HPDSolve");
#endif
    Cholesky( uplo, A );
    cholesky::SolveAfter( uplo, orientation, A, B );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
HPDSolve
( UpperOrLower uplo, Orientation orientation, 
  DistMatrix<F>& A, DistMatrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("HPDSolve");
#endif
    Cholesky( uplo, A );
    cholesky::SolveAfter( uplo, orientation, A, B );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef LAPACK_HPDSOLVE_HPP
